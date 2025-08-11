#!/usr/bin/env python3
"""
Biomedical Disease Tagger and Insight Tool
- Tries disease_map.json and disease_map2.json
- Validates JSON structure
- Caches WHO ICD-11 lookups per disease name
Requirements:
    pip install streamlit pandas requests
"""
import streamlit as st
import pandas as pd
import json
import logging
import requests
from typing import Dict, List, Tuple, Any
from pathlib import Path

# Logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

WHO_API_BASE = "https://id.who.int/icd/entity/search"

@st.cache_data
def load_disease_map() -> Dict[str, Dict[str, Any]]:
    """Try to load disease map from common filenames and normalize it to dict[str, dict]."""
    filenames = ["disease_map.json", "disease_map2.json"]
    for fname in filenames:
        p = Path(fname)
        if p.exists():
            try:
                with open(p, "r", encoding="utf-8") as f:
                    raw = json.load(f)
                if isinstance(raw, dict):
                    clean: Dict[str, Dict[str, Any]] = {}
                    for k, v in raw.items():
                        if not isinstance(v, dict):
                            logger.warning("Skipping entry %r because it's not a dict", k)
                            continue
                        clean[k.lower()] = v
                    logger.info("Loaded %s with %d valid entries", fname, len(clean))
                    return clean
                else:
                    logger.warning("%s exists but JSON root is not an object/dict", fname)
            except Exception as e:
                logger.error("Error loading %s: %s", fname, e)
    logger.warning("No disease_map file found or valid – returning empty mapping.")
    return {}

@st.cache_data
def fetch_icd11_code_from_who(disease_name: str) -> str:
    """Fetch ICD-11 code from WHO API. Cached per disease_name to avoid repeated calls."""
    api_key = st.secrets.get("WHO_API_KEY") if hasattr(st, "secrets") else None
    if not api_key:
        logger.warning("WHO_API_KEY not set in streamlit secrets. Skipping WHO lookup.")
        return "API key missing"

    try:
        headers = {"API-Version": "v2", "Authorization": f"Bearer {api_key}"}
        params = {"q": disease_name}
        r = requests.get(WHO_API_BASE, headers=headers, params=params, timeout=10)
        r.raise_for_status()
        data = r.json()
        if isinstance(data, dict) and "destinationEntities" in data and data["destinationEntities"]:
            first = data["destinationEntities"][0]
            entity_id = first.get("id", "")
            code = entity_id.split("/")[-1] if entity_id else "Not found"
            logger.info("WHO lookup %s -> %s", disease_name, code)
            return code
        return "Not found"
    except Exception as e:
        logger.error("Error fetching ICD-11 for %s: %s", disease_name, e)
        return "Error"

class DiseaseTagger:
    def __init__(self, disease_map: Dict[str, dict] = None):
        self.disease_map = disease_map or {}

    def tag_in_text(self, text: str) -> List[Tuple[str, str, str]]:
        text_lower = text.lower()
        matches: List[Tuple[str, str, str]] = []
        for disease, info in self.disease_map.items():
            if not isinstance(info, dict):
                continue
            if disease in text_lower:
                icd_code = info.get("icd11_code")
                if not icd_code or icd_code == "PLACEHOLDER":
                    icd_code = fetch_icd11_code_from_who(disease)
                matches.append((disease, info.get("therapeutic_area", "Not available"), icd_code))
        return matches

def create_streamlit_app():
    st.set_page_config(page_title="🧠 Disease Tagger with ICD-11", page_icon="🧠", layout="wide")

    # Custom CSS
    st.markdown("""
    <style>
        .main { background-color: #f7f9fc; }
        .stButton>button {
            background-color: #4CAF50; color: white; border-radius: 8px;
            padding: 0.5em 1em; font-weight: bold;
        }
        .stButton>button:hover { background-color: #45a049; }
        .section {
            background-color: white; padding: 20px; border-radius: 12px;
            box-shadow: 0 2px 6px rgba(0,0,0,0.1); margin-bottom: 20px;
        }
        h1 { text-align: center; color: #2C3E50; }
    </style>
    """, unsafe_allow_html=True)

    st.markdown("<h1>🧠 Disease Tagger with ICD-11</h1>", unsafe_allow_html=True)
    st.markdown("<p style='text-align:center;color:gray;'>Analyze biomedical text and view standardized ICD-11 codes (auto-updated via WHO API)</p>", unsafe_allow_html=True)

    # Load data
    disease_map = load_disease_map()
    tagger = DiseaseTagger(disease_map)

    # Sidebar
    with st.sidebar:
        st.header("Developer Tools")
        st.button("🔄 Refresh WHO ICD Cache")

    # Disease table
    st.markdown("<div class='section'>", unsafe_allow_html=True)
    st.subheader("📋 Available Diseases")
    if not disease_map:
        st.warning("⚠ No valid disease_map found. Place `disease_map.json` or `disease_map2.json` in the app folder.")
    else:
        rows = []
        for disease, info in disease_map.items():
            icd_code = info.get("icd11_code")
            if not icd_code or icd_code == "PLACEHOLDER":
                icd_code = fetch_icd11_code_from_who(disease)
            rows.append({
                "Disease": disease.title(),
                "Therapeutic Area": info.get("therapeutic_area", "Not available"),
                "ICD-11 Code": icd_code
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True)
    st.markdown("</div>", unsafe_allow_html=True)

    # Text tagging
    st.markdown("<div class='section'>", unsafe_allow_html=True)
    st.subheader("📝 Tag Your Text")
    user_text = st.text_area("Paste biomedical text here:", height=150, placeholder="e.g., Patient diagnosed with diabetes and hypertension.")
    if user_text:
        matches = tagger.tag_in_text(user_text)
        if matches:
            st.success(f"✅ Detected {len(matches)} disease(s):")
            for disease, area, icd in matches:
                st.write(f"**{disease.title()}** — {area} (ICD-11: `{icd}`)")
        else:
            st.warning("No known diseases detected.")
    st.markdown("</div>", unsafe_allow_html=True)

def main():
    create_streamlit_app()

if __name__ == "__main__":
    main()
