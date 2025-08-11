#!/usr/bin/env python3
"""
Biomedical Disease Tagger and Insight Tool
==========================================

This script provides functionality to:
1. Tag diseases in biomedical text
2. Map diseases to therapeutic areas
3. Fetch related PubMed abstracts
4. Provide a Streamlit web interface with ICD-11 codes

Requirements:
pip install streamlit pandas biopython spacy
python -m spacy download en_core_web_sm

Usage:
streamlit run disease_tagger.py
"""

import streamlit as st
import pandas as pd
import json
import logging
from difflib import get_close_matches
from typing import Dict, List, Tuple

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Streamlit cache functions
@st.cache_data
def load_disease_map():
    """Load disease mapping with ICD-11 codes from JSON."""
    try:
        with open("disease_map.json", "r") as f:
            return json.load(f)
    except Exception as e:
        logger.error(f"Error loading disease_map.json: {e}")
        return {}

class DiseaseTagger:
    def __init__(self, disease_map: Dict[str, dict] = None):
        self.disease_map = disease_map or load_disease_map()

    def tag_in_text(self, text: str) -> List[Tuple[str, str, str]]:
        text_lower = text.lower()
        matches = []
        for disease, info in self.disease_map.items():
            if disease in text_lower:
                matches.append((disease, info.get("therapeutic_area", "Not available"), info.get("icd11_code", "Not available")))
        return matches

def create_streamlit_app():
    st.set_page_config(page_title="Disease Tagger + ICD-11", layout="wide")
    st.title("🧠 Disease Tagger & ICD-11 Enhanced")
    st.markdown("Analyze biomedical text and view standardized ICD-11 codes.")

    # Load disease map
    disease_map = load_disease_map()
    tagger = DiseaseTagger(disease_map)

    # Display full disease table
    st.subheader("Available Diseases")
    rows = [
        {
            "Disease": disease.title(),
            "Therapeutic Area": info["therapeutic_area"],
            "ICD-11 Code": info["icd11_code"]
        }
        for disease, info in disease_map.items()
    ]
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
