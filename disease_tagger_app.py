#!/usr/bin/env python3
"""
Biomedical Disease Tagger and Insight Tool
==========================================

This script provides functionality to:
1. Tag diseases in biomedical text
2. Map diseases to therapeutic areas
3. Fetch related PubMed abstracts
4. Provide a Streamlit web interface

Requirements:
pip install streamlit pandas biopython spacy
python -m spacy download en_core_web_sm

Usage:
python disease_tagger.py
or
streamlit run disease_tagger.py
"""

import streamlit as st
import pandas as pd
import json
import os
import re
from typing import Dict, List, Tuple, Optional
from difflib import get_close_matches
import logging
import sys

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# ==============================
# Disease to Therapeutic Area Mapping
# ==============================

DISEASE_AREA_MAP = {
    "diabetes": "Endocrinology",
    "hypothyroidism": "Endocrinology",
    "cushing's syndrome": "Endocrinology",
    "acromegaly": "Endocrinology",

    "breast cancer": "Oncology",
    "lung cancer": "Oncology",
    "prostate cancer": "Oncology",
    "leukemia": "Oncology",
    "lymphoma": "Oncology",
    "glioblastoma": "Oncology",

    "asthma": "Pulmonology",
    "chronic obstructive pulmonary disease": "Pulmonology",
    "tuberculosis": "Pulmonology",
    "pulmonary fibrosis": "Pulmonology",

    "epilepsy": "Neurology",
    "alzheimer's disease": "Neurology",
    "parkinson's disease": "Neurology",
    "multiple sclerosis": "Neurology",
    "migraine": "Neurology",

    "hypertension": "Cardiology",
    "myocardial infarction": "Cardiology",
    "arrhythmia": "Cardiology",
    "congestive heart failure": "Cardiology",
    "stroke": "Cardiology",

    "crohn's disease": "Gastroenterology",
    "ulcerative colitis": "Gastroenterology",
    "hepatitis b": "Gastroenterology",
    "irritable bowel syndrome": "Gastroenterology",

    "covid-19": "Infectious Disease",
    "malaria": "Infectious Disease",
    "dengue": "Infectious Disease",
    "hiv": "Infectious Disease",

    "psoriasis": "Dermatology",
    "eczema": "Dermatology",
    "acne": "Dermatology",

    "chronic kidney disease": "Nephrology",
    "glomerulonephritis": "Nephrology",

    "rheumatoid arthritis": "Rheumatology",
    "lupus": "Rheumatology",

    "depression": "Psychiatry",
    "schizophrenia": "Psychiatry",
    "anxiety": "Psychiatry",
    "bipolar disorder": "Psychiatry"
}

# ==============================
# Sample Sentences
# ==============================

SAMPLE_SENTENCES = [
    "This study explores treatment options for diabetes in adults.",
    "Breast cancer patients were given a new monoclonal antibody.",
    "Epilepsy is a chronic condition affecting the brain.",
    "The prevalence of hypertension has increased among the elderly.",
    "COVID-19 vaccines showed high efficacy in preventing severe infection."
]

# ==============================
# Caching
# ==============================

@st.cache_data
def load_disease_map():
    return DISEASE_AREA_MAP


@st.cache_data
def save_disease_map_once():
    try:
        with open("disease_map.json", "w") as f:
            json.dump(DISEASE_AREA_MAP, f, indent=2)
        return True
    except Exception as e:
        logger.error(f"Error saving disease map: {e}")
        return False


# ==============================
# Disease Tagger Class
# ==============================

class DiseaseTagger:

    def __init__(self, disease_map: Dict[str, str] = None):
        self.disease_map = disease_map or load_disease_map()

    def tag_diseases_in_text(self, text: str) -> List[Tuple[str, str]]:
        found = []
        text_lower = text.lower()
        for disease in self.disease_map:
            if disease in text_lower:
                found.append((disease, self.disease_map[disease]))
        return found

    def find_similar_diseases(self, query: str, n: int = 3, cutoff: float = 0.4):
        return get_close_matches(query.lower(), self.disease_map.keys(), n=n, cutoff=cutoff)

    def process_sentences(self, sentences: List[str]) -> pd.DataFrame:
        results = []
        for text in sentences:
            diseases = self.tag_diseases_in_text(text)
            if diseases:
                for disease, area in diseases:
                    results.append({
                        "text": text,
                        "disease": disease,
                        "therapeutic_area": area
                    })
            else:
                results.append({
                    "text": text,
                    "disease": "Not found",
                    "therapeutic_area": "Unknown"
                })
        return pd.DataFrame(results)


# ==============================
# PubMed Fetcher
# ==============================

class PubMedFetcher:

    def __init__(self, email="your_email@example.com"):
        try:
            from Bio import Entrez
            Entrez.email = email
            self.entrez = Entrez
            self.available = True
        except:
            self.available = False

    def fetch_abstracts(self, query: str, max_results=3):
        if not self.available:
            return "Biopython not installed."

        try:
            handle = self.entrez.esearch(db="pubmed", term=query, retmax=max_results)
            record = self.entrez.read(handle)
            ids = ",".join(record["IdList"])
            fetch = self.entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="text")
            return fetch.read()
        except Exception as e:
            return str(e)


# ==============================
# CLI DEMO
# ==============================

def run_command_line_demo():
    print("🧠 Disease Tagger CLI Demo")
    tagger = DiseaseTagger()
    pubmed = PubMedFetcher()

    while True:
        text = input("Enter sentence (or quit): ")
        if text.lower() in ["quit", "exit"]:
            break

        diseases = tagger.tag_diseases_in_text(text)

        if diseases:
            for disease, area in diseases:
                print(f"Detected: {disease} → {area}")
                abstracts = pubmed.fetch_abstracts(disease)
                print(abstracts[:500])
        else:
            print("No disease found.")


# ==============================
# STREAMLIT APP
# ==============================

def create_streamlit_app():

    st.set_page_config(
        page_title="Disease Tagger & Insight App",
        page_icon="🧠",
        layout="wide"
    )

    st.title("🧠 Disease Tagger & Insight App")
    st.markdown("Analyze biomedical text and discover therapeutic insights")

    tagger = DiseaseTagger()
    pubmed = PubMedFetcher()

    text = st.text_area("Enter biomedical text:")

    if text:
        diseases = tagger.tag_diseases_in_text(text)

        if diseases:
            for i, (disease, area) in enumerate(diseases):
                st.success(f"{disease.title()} → {area}")
                if st.button("Fetch Abstracts", key=i):
                    abstracts = pubmed.fetch_abstracts(disease)
                    st.text_area("PubMed Results", abstracts, height=300)
        else:
            st.warning("No disease detected.")


# ==============================
# MAIN EXECUTION (UNCHANGED LOGIC)
# ==============================

def main():
    save_disease_map_once()

    if len(sys.argv) > 1 and sys.argv[1] == "--cli":
        run_command_line_demo()
    else:
        create_streamlit_app()


if __name__ == "__main__":
    main()
