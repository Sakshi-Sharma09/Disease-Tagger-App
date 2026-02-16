#!/usr/bin/env python3
"""
Biomedical Disease Tagger and Insight Tool
"""

import streamlit as st
import pandas as pd
import json
import os
import re
from typing import Dict, List, Tuple
from difflib import get_close_matches
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# ==============================
# DISEASE MAP
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
# CACHE
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
    except Exception:
        return False

# ==============================
# CORE CLASSES
# ==============================

class DiseaseTagger:
    def __init__(self):
        self.disease_map = load_disease_map()

    def tag_diseases_in_text(self, text: str):
        results = []
        text_lower = text.lower()
        for disease in self.disease_map:
            if disease in text_lower:
                results.append((disease, self.disease_map[disease]))
        return results

    def find_similar_diseases(self, query: str, n=3, cutoff=0.4):
        return get_close_matches(query.lower(), self.disease_map.keys(), n=n, cutoff=cutoff)

    def process_sentences(self, sentences: List[str]):
        data = []
        for text in sentences:
            diseases = self.tag_diseases_in_text(text)
            if diseases:
                for disease, area in diseases:
                    data.append({
                        "text": text,
                        "disease": disease,
                        "therapeutic_area": area
                    })
            else:
                data.append({
                    "text": text,
                    "disease": "Not found",
                    "therapeutic_area": "Unknown"
                })
        return pd.DataFrame(data)

class PubMedFetcher:
    def __init__(self, email="your_email@example.com"):
        self.available = False
        try:
            from Bio import Entrez
            Entrez.email = email
            self.entrez = Entrez
            self.available = True
        except:
            pass

    def fetch_abstracts(self, query, max_results=3):
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
# STREAMLIT APP
# ==============================

def create_streamlit_app():

    st.set_page_config(
        page_title="Disease Tagger & Insight App",
        page_icon="🧠",
        layout="wide"
    )

    if "logged_in" not in st.session_state:
        st.session_state.logged_in = False

    if not st.session_state.logged_in:

        st.title("🧠 Disease Tagger & Insight App")
        st.subheader("🔐 Secure Access Portal")

        if os.path.exists("users.json"):
            with open("users.json", "r") as f:
                USERS = json.load(f)
        else:
            USERS = {}

        option = st.selectbox("Select Option", ["Login", "Sign Up"])

        username = st.text_input("Username")
        password = st.text_input("Password", type="password")

        if option == "Login":
            if st.button("Login"):
                if username in USERS and USERS[username] == password:
                    st.session_state.logged_in = True
                    st.session_state.username = username
                    st.rerun()
                else:
                    st.error("Invalid credentials")

        if option == "Sign Up":
            if st.button("Create Account"):
                USERS[username] = password
                with open("users.json", "w") as f:
                    json.dump(USERS, f)
                st.success("Account created!")

        st.stop()

    # AFTER LOGIN

    with st.sidebar:
        st.success(f"Logged in as {st.session_state.username}")
        if st.button("Logout"):
            st.session_state.logged_in = False
            st.rerun()

    st.title("🧠 Disease Tagger & Insight App")

    tagger = DiseaseTagger()
    pubmed = PubMedFetcher()

    text = st.text_area("Enter biomedical text")

    if text:
        diseases = tagger.tag_diseases_in_text(text)

        if diseases:
            for disease, area in diseases:
                st.success(f"{disease.title()} → {area}")
                if st.button(f"Fetch Abstracts for {disease}"):
                    abstracts = pubmed.fetch_abstracts(disease)
                    st.text_area("PubMed Results", abstracts, height=300)
        else:
            st.warning("No disease detected.")

# ==============================
# RUN STREAMLIT DIRECTLY
# ==============================

save_disease_map_once()
create_streamlit_app()
