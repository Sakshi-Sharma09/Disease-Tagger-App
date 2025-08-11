#!/usr/bin/env python3
"""
ICD-11 Disease Tagger — async + batch + UI
Requirements:
    pip install streamlit pandas aiohttp
Usage:
    streamlit run ICD-11.py
Notes:
 - Put disease_map.json or disease_map2.json next to this file.
 - Add WHO API token to Streamlit secrets as WHO_API_KEY for live lookups.
"""
import asyncio
import aiohttp
import streamlit as st
import pandas as pd
import json
import logging
import re
from pathlib import Path
from typing import Dict, List, Tuple, Any
from datetime import datetime

# ---------- Config ----------
DISEASE_FILES = ["disease_map.json", "disease_map2.json"]
ICD_CACHE_FILE = "icd11_cache.json"
WHO_SEARCH_ENDPOINT = "https://id.who.int/icd/entity/search"  # simple search endpoint
MAX_CONCURRENT_REQUESTS = 8  # tune for concurrency
# ----------------------------

# Logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("icd11-tagger")

# ---------- Helpers: load/save ----------
@st.cache_data
def load_disease_map() -> Dict[str, Dict[str, Any]]:
    """Load disease mapping from one of expected files and normalize keys to lowercase."""
    for fname in DISEASE_FILES:
        p = Path(fname)
        if p.exists():
            try:
                raw = json.loads(p.read_text(encoding="utf-8"))
                if isinstance(raw, dict):
                    clean: Dict[str, Dict[str, Any]] = {}
                    for k, v in raw.items():
                        if isinstance(v, dict):
                            clean[k.lower()] = v
                        else:
                            logger.warning("Skipping entry %r because it's not an object/dict", k)
                    return clean
                else:
                    logger.warning("%s exists but root is not a dict", fname)
            except Exception as e:
                logger.error("Error reading %s: %s", fname, e)
    return {}

def read_icd_cache() -> Dict[str, str]:
    p = Path(ICD_CACHE_FILE)
    if not p.exists():
        return {}
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return {}

def write_icd_cache(cache: Dict[str, str]) -> None:
    try:
        Path(ICD_CACHE_FILE).write_text(json.dumps(cache, indent=2), encoding="utf-8")
    except Exception as e:
        logger.error("Failed to write ICD cache: %s", e)

# ---------- WHO lookup (async) ----------
async def async_fetch_single(session: aiohttp.ClientSession, disease: str, api_key: str) -> Tuple[str, str]:
    """Perform one WHO query; return (disease, code)."""
    # The WHO API accepts q param on some endpoints; using search endpoint here.
    params = {"q": disease}
    headers = {"API-Version": "v2"}
    if api_key:
        headers["Authorization"] = f"Bearer {api_key}"
    try:
        async with session.get(WHO_SEARCH_ENDPOINT, params=params, headers=headers, timeout=15) as resp:
            if resp.status != 200:
                logger.debug("WHO lookup %s returned %s", disease, resp.status)
                return disease, "Not found"
            data = await resp.json()
            # destinationEntities is the common key; try multiple fallbacks defensively
            dest = data.get("destinationEntities") or data.get("destinationEntities", [])
            if dest:
                first = dest[0]
                entity_id = first.get("id") or first.get("theCode") or ""
                code = entity_id.split("/")[-1] if entity_id else "Not found"
                return disease, code
    except Exception as e:
        logger.debug("WHO lookup failed for %s: %s", disease, e)
    return disease, "Error"

async def async_fetch_all(diseases: List[str], api_key: str, progress_callback=None) -> Dict[str, str]:
    """Fetch codes concurrently for all diseases; returns mapping disease->code."""
    cache = read_icd_cache()
    # only query those not in cache
    to_query = [d for d in diseases if d not in cache]
    results: Dict[str, str] = dict(cache)  # start with existing cache

    if not to_query:
        return results

    conn = aiohttp.TCPConnector(limit_per_host=MAX_CONCURRENT_REQUESTS)
    timeout = aiohttp.ClientTimeout(total=60)
    semaphore = asyncio.Semaphore(MAX_CONCURRENT_REQUESTS)

    async def sem_fetch(session, d):
        async with semaphore:
            return await async_fetch_single(session, d, api_key)

    async with aiohttp.ClientSession(connector=conn, timeout=timeout) as session:
        tasks = [asyncio.create_task(sem_fetch(session, d)) for d in to_query]
        for i, t in enumerate(asyncio.as_completed(tasks), start=1):
            pair = await t
            disease, code = pair
            results[disease] = code
            # update persisted cache incrementally for resilience
            try:
                write_icd_cache(results)
            except Exception:
                pass
            if progress_callback:
                progress_callback(len(results), len(diseases))
    return results

# ---------- Tagger ----------
class DiseaseTagger:
    def __init__(self, disease_map: Dict[str, dict]):
        self.disease_map = disease_map or {}

    def tag_in_text(self, text: str) -> List[Tuple[str, str, str]]:
        text_lower = text.lower()
        matches: List[Tuple[str, str, str]] = []
        for disease, info in self.disease_map.items():
            if not isinstance(info, dict):
                continue
            # word boundary match for safety
            if re.search(rf"\b{re.escape(disease)}\b", text_lower):
                icd_code = info.get("icd11_code")
                matches.append((disease, info.get("therapeutic_area", "Not available"), icd_code or "PLACEHOLDER"))
        return matches

# ---------- UI ----------
def create_streamlit_app():
    st.set_page_config(page_title="🧠 Disease Tagger with ICD-11", page_icon="🧠", layout="wide")
    # ----- styles -----
    st.markdown(
        """
    <style>
      .main { background: linear-gradient(180deg,#f7fbff 0%, #ffffff 100%); padding: 1rem 2rem;}
      .card { background: white; padding: 18px; border-radius: 12px; box-shadow: 0 6px 18px rgba(12,24,48,0.06); }
      .big-metric { font-size: 1.3rem; font-weight:700; color:#0f172a; }
      .muted { color:#6b7280; }
      .section { margin-bottom: 20px; }
      .small { font-size:0.85rem; color:#6b7280; }
    </style>
    """, unsafe_allow_html=True
    )

    st.markdown("<h1 style='text-align:center;'>🧠 Disease Tagger with ICD-11</h1>", unsafe_allow_html=True)
    st.markdown("<p style='text-align:center;color:gray'>Tag biomedical text, batch-process files and fetch ICD-11 codes asynchronously</p>", unsafe_allow_html=True)

    # load mappings and cache
    disease_map = load_disease_map()
    tagger = DiseaseTagger(disease_map)
    icd_cache = read_icd_cache()
    who_api_key = st.secrets.get("WHO_API_KEY") if hasattr(st, "secrets") else None

    # Sidebar tools
    with st.sidebar:
        st.header("⚙️ Developer Tools")
        if st.button("🔄 Clear ICD Cache"):
            write_icd_cache({})
            st.experimental_rerun()
        st.markdown("---")
        st.markdown("**WHO API Key**")
        if who_api_key:
            st.success("WHO API key loaded (from secrets).")
        else:
            st.warning("No WHO API key in `st.secrets`. Live lookups disabled.")
        st.caption("Add `WHO_API_KEY` in Streamlit Cloud → Manage App → Secrets")

    # Top stats cards
    col1, col2, col3, col4 = st.columns([1.2,1,1,1])
    total_diseases = len(disease_map)
    cached_codes = len(icd_cache)
    col1.markdown("<div class='card'><div class='big-metric'>📚 Diseases</div><div class='muted'>In DB</div><div style='font-size:20px;font-weight:700'>{}</div></div>".format(total_diseases), unsafe_allow_html=True)
    col2.markdown("<div class='card'><div class='big-metric'>🔁 Cached</div><div class='muted'>ICD-11 codes</div><div style='font-size:20px;font-weight:700'>{}</div></div>".format(cached_codes), unsafe_allow_html=True)
    last_mod = Path(ICD_CACHE_FILE).stat().st_mtime if Path(ICD_CACHE_FILE).exists() else None
    last_mod_str = datetime.fromtimestamp(last_mod).strftime("%Y-%m-%d %H:%M:%S") if last_mod else "—"
    col3.markdown("<div class='card'><div class='big-metric'>🕘 Last cache</div><div class='muted'>Updated</div><div style='font-size:18px'>{}</div></div>".format(last_mod_str), unsafe_allow_html=True)
    col4.markdown("<div class='card'><div class='big-metric'>⚙️ WHO</div><div class='muted'>Key present</div><div style='font-size:18px'>{}</div></div>".format("Yes" if who_api_key else "No"), unsafe_allow_html=True)

    st.markdown("---")

    # Tabs: Single / Batch / Table / Manage
    tab1, tab2, tab3, tab4 = st.tabs(["🔍 Single Tag", "📂 Batch Process", "📋 Disease Table", "🛠 Manage"])

    # ---- Tab 1: Single Tag ----
    with tab1:
        st.markdown("<div class='card section'>", unsafe_allow_html=True)
        st.subheader("Paste text to tag diseases")
        single_text = st.text_area("Enter biomedical text...", height=150, placeholder="e.g., Patient presents with diabetes and hypertension.")
        if st.button("🔎 Tag Text"):
            if not single_text.strip():
                st.warning("Please enter text to analyze.")
            else:
                matches = tagger.tag_in_text(single_text)
                if not matches:
                    st.info("No known diseases found in the input text.")
                else:
                    # ensure ICD codes loaded for PLACEHOLDER entries (may be async)
                    output_rows = []
                    for d, area, code in matches:
                        if (not code or code == "PLACEHOLDER") and who_api_key:
                            # attempt quick synchronous check in cache
                            cached = icd_cache.get(d)
                            code = cached or "Fetching..."
                        output_rows.append({"Disease": d.title(), "Therapeutic Area": area, "ICD-11": code})
                    df_out = pd.DataFrame(output_rows)
                    st.dataframe(df_out, use_container_width=True)
                    # Offer to fetch missing codes async
                    missing = [r["Disease"].lower() for r in output_rows if r["ICD-11"] in ("PLACEHOLDER", "Fetching...")]
                    if missing and who_api_key:
                        if st.button("⚡ Fetch missing ICD-11 codes (async)"):
                            progress = st.progress(0)
                            status_text = st.empty()
                            async def runner():
                                def cb(done, total):
                                    progress.progress(int(done/total*100))
                                    status_text.info(f"Fetched {done}/{total}")
                                res = await async_fetch_all(missing, who_api_key, progress_callback=lambda a,b: cb(a,b))
                                return res
                            res_map = asyncio.run(runner())
                            # show results
                            rows2 = []
                            for d, area, _ in matches:
                                code = res_map.get(d) or icd_cache.get(d) or "Not found"
                                rows2.append({"Disease": d.title(), "Therapeutic Area": area, "ICD-11": code})
                            st.success("✅ Done")
                            st.dataframe(pd.DataFrame(rows2), use_container_width=True)
        st.markdown("</div>", unsafe_allow_html=True)

    # ---- Tab 2: Batch ----
    with tab2:
        st.markdown("<div class='card section'>", unsafe_allow_html=True)
        st.subheader("Batch process TXT or CSV (column-wise)")
        uploaded = st.file_uploader("Upload TXT or CSV", type=["txt","csv"])
        concurrency = st.slider("Concurrency (parallel WHO requests)", min_value=2, max_value=20, value=8)
        if uploaded:
            content = uploaded.read()
            if uploaded.name.lower().endswith(".txt"):
                try:
                    text = content.decode("utf-8")
                except Exception:
                    text = content.decode("latin-1")
                # simple whole-file tag
                matches = tagger.tag_in_text(text)
                out_df = pd.DataFrame(matches, columns=["disease","therapeutic_area","icd11_code"])
                st.write(f"Found {len(out_df)} matches in TXT")
                st.dataframe(out_df, use_container_width=True)
            else:
                # csv
                try:
                    df_in = pd.read_csv(uploaded)
                except Exception:
                    st.error("Unable to read CSV. Ensure it's valid.")
                    df_in = None
                if df_in is not None:
                    st.write("Columns detected:", list(df_in.columns))
                    # flatten all text fields into single list to extract candidate diseases
                    text_cells = []
                    for _, row in df_in.iterrows():
                        # join row fields
                        row_text = " ".join([str(x) for x in row.values if pd.notna(x)])
                        text_cells.append(row_text)
                    # find matches per row
                    progress_bar = st.progress(0)
                    status = st.empty()
                    all_diseases = list(disease_map.keys())
                    results_rows = []
                    # We'll find candidate diseases first, then fetch missing ICDs async
                    for idx, row_text in enumerate(text_cells):
                        found = tagger.tag_in_text(row_text)
                        for d, area, code in found:
                            results_rows.append({"row": idx, "Disease": d.title(), "Therapeutic Area": area, "ICD-11": code or "PLACEHOLDER"})
                        progress_bar.progress(int((idx+1)/len(text_cells)*100))
                    st.write(f"Pre-scan found {len(results_rows)} disease mentions across file.")
                    df_res = pd.DataFrame(results_rows)
                    if df_res.empty:
                        st.info("No disease mentions found in file.")
                    else:
                        st.dataframe(df_res, use_container_width=True)
                        # ask to run async WHO fetch for missing codes
                        missing_list = list({r["Disease"].lower() for r in results_rows if r["ICD-11"] in (None, "", "PLACEHOLDER")})
                        if missing_list and who_api_key:
                            if st.button(f"⚡ Fetch ICD-11 for {len(missing_list)} missing items (async)"):
                                # run async fetch with user-controlled concurrency
                                global MAX_CONCURRENT_REQUESTS
                                MAX_CONCURRENT_REQUESTS = concurrency
                                prog = st.progress(0)
                                info = st.empty()
                                def cb(done, total):
                                    prog.progress(int(done/total*100))
                                    info.info(f"Fetched {done}/{total}")
                                async def run_batch():
                                    return await async_fetch_all(missing_list, who_api_key, progress_callback=lambda a,b: cb(a,b))
                                res_map = asyncio.run(run_batch())
                                # merge results into df_res
                                df_res["ICD-11"] = df_res.apply(lambda r: res_map.get(r["Disease"].lower(), r["ICD-11"]), axis=1)
                                st.success("Batch ICD fetch done.")
                                st.dataframe(df_res, use_container_width=True)
                                # offer CSV download
                                st.download_button("⬇️ Download results CSV", df_res.to_csv(index=False).encode("utf-8"), file_name="batch_icd_results.csv")
        st.markdown("</div>", unsafe_allow_html=True)

    # ---- Tab 3: Table ----
    with tab3:
        st.markdown("<div class='card section'>", unsafe_allow_html=True)
        st.subheader("Full disease table view")
        if not disease_map:
            st.warning("No disease_map found. Add disease_map.json or disease_map2.json to app folder.")
        else:
            # prepare table with cached codes where missing
            combined_cache = read_icd_cache()
            rows = []
            for disease, info in disease_map.items():
                code = info.get("icd11_code") or combined_cache.get(disease) or ""
                rows.append({"Disease": disease.title(), "Therapeutic Area": info.get("therapeutic_area","Not available"), "ICD-11": code})
            df_table = pd.DataFrame(rows)
            st.dataframe(df_table, use_container_width=True)
            if st.button("🔁 Refresh missing ICD (async, cached)"):
                # fetch all missing
                missing = [r["Disease"].lower() for r in rows if not r["ICD-11"]]
                if not missing:
                    st.info("No missing codes to fetch.")
                elif not who_api_key:
                    st.warning("WHO API key missing — put WHO_API_KEY into Streamlit secrets to enable.")
                else:
                    p = st.progress(0)
                    info = st.empty()
                    def cb(d,t): p.progress(int(d/t*100)); info.info(f"{d}/{t}")
                    res_map = asyncio.run(async_fetch_all(missing, who_api_key, progress_callback=lambda a,b: cb(a,b)))
                    st.success("Done — cache updated.")
                    st.experimental_rerun()
        st.markdown("</div>", unsafe_allow_html=True)

    # ---- Tab 4: Manage ----
    with tab4:
        st.markdown("<div class='card section'>", unsafe_allow_html=True)
        st.subheader("Manage & Export")
        st.markdown("You can persist newly fetched ICD-11 codes back to `disease_map.json` (will overwrite file).")
        if st.checkbox("Show sample disease_map.json structure"):
            sample = {
                "diabetes": {"therapeutic_area":"Endocrinology","icd11_code":"5A11"},
                "asthma": {"therapeutic_area":"Pulmonology","icd11_code":"CA23"}
            }
            st.json(sample)
        if st.button("🔁 Save cache back to disease_map.json"):
            # read both structures, merge cache into disease_map and write to disease_map.json
            cache = read_icd_cache()
            if not cache:
                st.warning("Cache empty — nothing to write.")
            else:
                # load existing (prefer first file in DISEASE_FILES or create new)
                target = Path(DISEASE_FILES[0])
                existing = {}
                if target.exists():
                    try:
                        existing = json.loads(target.read_text(encoding="utf-8"))
                    except Exception as e:
                        st.error(f"Failed reading {target}: {e}")
                        existing = {}
                # update
                for disease, code in cache.items():
                    key = disease.lower()
                    entry = existing.get(key, {})
                    if not isinstance(entry, dict):
                        entry = {}
                    entry["icd11_code"] = code
                    existing[key] = entry
                try:
                    target.write_text(json.dumps(existing, indent=2), encoding="utf-8")
                    st.success(f"Saved {len(cache)} ICD entries to {target.name}")
                except Exception as e:
                    st.error(f"Failed to write file: {e}")
        st.markdown("</div>", unsafe_allow_html=True)

# ---------- Main ----------
if __name__ == "__main__":
    create_streamlit_app()
