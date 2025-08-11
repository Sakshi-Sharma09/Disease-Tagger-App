#!/usr/bin/env python3
"""
Biomedical Disease Tagger and Insight Tool
==========================================

Now includes automatic ICD-11 code fetching from WHO API if not found in local JSON.

Requirements:
pip install streamlit pandas requests biopython spacy
python -m spacy download en_core_web_sm
"""

import streamlit as st
import pandas as pd
import json
import logging
import requests
from difflib import get_close_matches
from typing import Dict, List, Tuple

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

WHO_API_BASE = "https://id.who.int/icd/entity/search"

@st.cache_data
def load_disease_map():
    """Load disease mapping with ICD-11 codes from JSON."""
    try:
        with open("disease_map.json", "r") as f:
            return json.load(f)
    except Exception as e:
        logger.error(f"Error loading disease_map.json: {e}")
        return {}

@st.cache_data
def fetch_icd11_code_from_who(disease_name: str) -> str:
    """Fetch ICD-11 code from WHO API."""
    api_key = st.secrets.get("WHO_API_KEY")
    if not api_key:
        logger.warning("WHO_API_KEY not found in secrets.")
        return "API key missing"

    try:
        headers = {"API-Version": "v2", "Authorization": f"Bearer {api_key}"}
        params = {"q": disease_name}
        r = requests.get(WHO_API_BASE, headers=headers, params=params, timeout=10)
        r.raise_for_status()
        data = r.json()

        if "destinationEntities" in data and data["destinationEntities"]:
            return data["destinationEntities"][0].get("id", "Not found").split("/")[-1]
        return "Not found"
    except Exception as e:
        logger.error(f"Error fetching ICD-11 for {disease_name}: {e}")
        return "Error"

class DiseaseTagger:
    def __init__(self, disease_map: Dict[str, dict] = None):
        self.disease_map = disease_map or load_disease_map()

    def tag_in_text(self, text: str) -> List[Tuple[str, str, str]]:
        text_lower = text.lower()
        matches = []
        for disease, info in self.disease_map.items():
            if disease in text_lower:
                icd_code = info.get("icd11_code")
                if not icd_code or icd_code == "PLACEHOLDER":
                    icd_code = fetch_icd11_code_from_who(disease)
                matches.append((disease, info.get("therapeutic_area", "Not available"), icd_code))
        return matches

def create_streamlit_app():
    st.set_page_config(page_title="Disease Tagger + ICD-11 Auto", layout="wide")
    st.title("🧠 Disease Tagger & ICD-11 Auto-Updater")
    st.markdown("Analyze biomedical text and view standardized ICD-11 codes (auto-updated via WHO API).")

    # Load disease map
    disease_map = load_disease_map()
    tagger = DiseaseTagger(disease_map)

    # Display full disease table
    st.subheader("Available Diseases")
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

    # Input and tagging
    st.subheader("Tag Your Text")
    user_text = st.text_area("Paste biomedical text here:", height=150)
    if user_text:
        matches = tagger.tag_in_text(user_text)
        if matches:
            st.success(f"Detected {len(matches)} disease(s):")
            for disease, area, icd in matches:
                st.write(f"**{disease.title()}** — {area} (ICD-11: {icd})")
        else:
            st.warning("No known diseases detected.")

def main():
    create_streamlit_app()

if __name__ == "__main__":
    main()
