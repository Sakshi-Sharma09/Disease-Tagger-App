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

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Disease to Therapeutic Area Mapping
DISEASE_AREA_MAP = {
    # Endocrinology
    "diabetes": "Endocrinology",
    "hypothyroidism": "Endocrinology",
    "cushing's syndrome": "Endocrinology",
    "acromegaly": "Endocrinology",

    # Oncology
    "breast cancer": "Oncology",
    "lung cancer": "Oncology",
    "prostate cancer": "Oncology",
    "leukemia": "Oncology",
    "lymphoma": "Oncology",
    "glioblastoma": "Oncology",

    # Pulmonology
    "asthma": "Pulmonology",
    "chronic obstructive pulmonary disease": "Pulmonology",
    "tuberculosis": "Pulmonology",
    "pulmonary fibrosis": "Pulmonology",

    # Neurology
    "epilepsy": "Neurology",
    "alzheimer's disease": "Neurology",
    "parkinson's disease": "Neurology",
    "multiple sclerosis": "Neurology",
    "migraine": "Neurology",

    # Cardiology
    "hypertension": "Cardiology",
    "myocardial infarction": "Cardiology",
    "arrhythmia": "Cardiology",
    "congestive heart failure": "Cardiology",
    "stroke": "Cardiology",

    # Gastroenterology
    "crohn's disease": "Gastroenterology",
    "ulcerative colitis": "Gastroenterology",
    "hepatitis b": "Gastroenterology",
    "irritable bowel syndrome": "Gastroenterology",

    # Infectious Diseases
    "covid-19": "Infectious Disease",
    "malaria": "Infectious Disease",
    "dengue": "Infectious Disease",
    "hiv": "Infectious Disease",

    # Dermatology
    "psoriasis": "Dermatology",
    "eczema": "Dermatology",
    "acne": "Dermatology",

    # Nephrology
    "chronic kidney disease": "Nephrology",
    "glomerulonephritis": "Nephrology",

    # Rheumatology
    "rheumatoid arthritis": "Rheumatology",
    "lupus": "Rheumatology",

    # Psychiatry
    "depression": "Psychiatry",
    "schizophrenia": "Psychiatry",
    "anxiety": "Psychiatry",
    "bipolar disorder": "Psychiatry"
}

# Sample biomedical sentences for testing
SAMPLE_SENTENCES = [
    "This study explores treatment options for diabetes in adults.",
    "Breast cancer patients were given a new monoclonal antibody.",
    "Epilepsy is a chronic condition affecting the brain.",
    "The prevalence of hypertension has increased among the elderly.",
    "COVID-19 vaccines showed high efficacy in preventing severe infection.",
    "Lung cancer is often diagnosed at a late stage, making treatment difficult.",
    "Migraine is a common neurological disorder with significant impact on quality of life.",
    "Crohn's disease can lead to inflammation throughout the gastrointestinal tract.",
    "Patients with psoriasis reported improvement after using the new biologic agent.",
    "Parkinson's disease is associated with tremors and motor dysfunction.",
    "Asthma exacerbations were reduced with inhaled corticosteroids.",
    "HIV remains a major global health issue, particularly in low-income countries.",
    "Chronic kidney disease requires early diagnosis and dietary modifications.",
    "Myocardial infarction is commonly known as a heart attack.",
    "Rheumatoid arthritis is an autoimmune disease that affects joints.",
    "Schizophrenia is characterized by hallucinations, delusions, and cognitive dysfunction.",
    "Alzheimer's disease leads to progressive memory loss in elderly patients.",
    "Patients with hepatitis B are at increased risk for liver cirrhosis.",
    "Eczema flares can be triggered by allergens and environmental factors.",
    "Dengue fever is transmitted by mosquitoes in tropical regions.",
]

# Streamlit caching functions
@st.cache_data
def load_disease_map():
    """Load and cache the disease map to prevent repeated execution."""
    return DISEASE_AREA_MAP

@st.cache_data
def save_disease_map_once():
    """Save disease map only once using Streamlit caching."""
    try:
        with open("disease_map.json", "w") as f:
            json.dump(DISEASE_AREA_MAP, f, indent=2)
        logger.info("Disease map saved to disease_map.json")
        return True
    except Exception as e:
        logger.error(f"Error saving disease map: {e}")
        return False

class DiseaseTagger:
    """Main class for disease tagging functionality."""
    
    def __init__(self, disease_map: Dict[str, str] = None):
        # Use cached disease map
        self.disease_map = disease_map or load_disease_map()
        
    def tag_diseases_in_text(self, text: str) -> List[Tuple[str, str]]:
        """
        Tag diseases found in the given text.
        
        Args:
            text (str): Input text to analyze
            
        Returns:
            List[Tuple[str, str]]: List of (disease, therapeutic_area) tuples
        """
        found_diseases = []
        text_lower = text.lower()
        
        for disease in self.disease_map:
            if disease in text_lower:
                therapeutic_area = self.disease_map[disease]
                found_diseases.append((disease, therapeutic_area))
                
        return found_diseases
    
    def find_similar_diseases(self, query: str, n: int = 3, cutoff: float = 0.4) -> List[str]:
        """
        Find diseases similar to the query using fuzzy matching.
        
        Args:
            query (str): Disease name to search for
            n (int): Number of matches to return
            cutoff (float): Similarity threshold
            
        Returns:
            List[str]: List of similar disease names
        """
        return get_close_matches(query.lower(), self.disease_map.keys(), n=n, cutoff=cutoff)
    
    def process_sentences(self, sentences: List[str]) -> pd.DataFrame:
        """
        Process multiple sentences and return results as DataFrame.
        
        Args:
            sentences (List[str]): List of sentences to process
            
        Returns:
            pd.DataFrame: DataFrame with text, disease, and therapeutic area columns
        """
        results = []
        
        for text in sentences:
            diseases = self.tag_diseases_in_text(text)
            
            if diseases:
                for disease, area in diseases:
                    results.append({
                        'text': text,
                        'disease': disease,
                        'therapeutic_area': area
                    })
            else:
                results.append({
                    'text': text,
                    'disease': 'Not found',
                    'therapeutic_area': 'Unknown'
                })
                
        return pd.DataFrame(results)

class PubMedFetcher:
    """Class for fetching PubMed abstracts."""
    
    def __init__(self, email: str = "your_email@example.com"):
        self.email = email
        try:
            from Bio import Entrez
            Entrez.email = self.email
            self.entrez = Entrez
            self.available = True
        except ImportError:
            logger.warning("Biopython not available. PubMed functionality disabled.")
            self.available = False
    
    def fetch_abstracts(self, query: str, max_results: int = 3) -> str:
        """
        Fetch PubMed abstracts for a given query.
        
        Args:
            query (str): Search query
            max_results (int): Maximum number of results to fetch
            
        Returns:
            str: Formatted abstracts or error message
        """
        if not self.available:
            return "PubMed functionality not available. Please install biopython: pip install biopython"
        
        try:
            handle = self.entrez.esearch(db="pubmed", term=query, retmax=max_results)
            record = self.entrez.read(handle)
            
            if not record["IdList"]:
                return f"No abstracts found for query: {query}"
            
            ids = ",".join(record["IdList"])
            fetch = self.entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="text")
            abstracts = fetch.read()
            
            return abstracts if abstracts else f"No abstract content found for query: {query}"
            
        except Exception as e:
            logger.error(f"Error fetching PubMed abstracts: {e}")
            return f"Error fetching abstracts: {e}"

def save_disease_map_to_file(filename: str = "disease_map.json"):
    """Save the disease map to a JSON file."""
    try:
        with open(filename, "w") as f:
            json.dump(DISEASE_AREA_MAP, f, indent=2)
        logger.info(f"Disease map saved to {filename}")
    except Exception as e:
        logger.error(f"Error saving disease map: {e}")

def load_disease_map_from_file(filename: str = "disease_map.json") -> Dict[str, str]:
    """Load disease map from a JSON file."""
    try:
        with open(filename, "r") as f:
            return json.load(f)
    except FileNotFoundError:
        logger.warning(f"File {filename} not found. Using default disease map.")
        return DISEASE_AREA_MAP
    except Exception as e:
        logger.error(f"Error loading disease map: {e}")
        return DISEASE_AREA_MAP

def run_command_line_demo():
    """Run a command-line demonstration of the disease tagger."""
    print("🧠 Disease Tagger - Command Line Demo")
    print("=" * 50)
    
    tagger = DiseaseTagger()
    pubmed = PubMedFetcher()
    
    # Process sample sentences
    print("\nProcessing sample sentences:")
    df = tagger.process_sentences(SAMPLE_SENTENCES)
    print(df.head(10))
    
    # Interactive mode
    while True:
        print("\n" + "=" * 50)
        user_input = input("Enter a biomedical sentence (or 'quit' to exit): ").strip()
        
        if user_input.lower() in ['quit', 'exit', 'q']:
            break
            
        if not user_input:
            continue
            
        diseases = tagger.tag_diseases_in_text(user_input)
        
        if diseases:
            for disease, area in diseases:
                print(f"\n🔬 Disease Detected: {disease.title()}")
                print(f"📋 Therapeutic Area: {area}")
                
                # Fetch PubMed abstracts
                print("📄 Fetching PubMed abstracts...")
                abstracts = pubmed.fetch_abstracts(disease, max_results=2)
                print(abstracts[:500] + "..." if len(abstracts) > 500 else abstracts)
        else:
            print("❌ No known diseases detected.")
            
            # Try fuzzy matching
            words = user_input.lower().split()
            for word in words:
                similar = tagger.find_similar_diseases(word)
                if similar:
                    print(f"💡 Did you mean: {', '.join(similar)}?")
                    break

def create_streamlit_app():
    """Create and run the Streamlit web application."""
    st.set_page_config(
        page_title="Disease Tagger & Insight Tool",
        page_icon="🧠",
        layout="wide"
    )
    # -------- LOGIN / SIGNUP SYSTEM --------
    if "logged_in" not in st.session_state:
        st.session_state.logged_in = False

    if not st.session_state.logged_in:

        st.markdown("<h1 style='text-align: center;'>🧠 Disease Tagger & Insight Tool</h1>", unsafe_allow_html=True)
        st.markdown("<h3 style='text-align: center;'>🔐 Secure Access Portal</h3>", unsafe_allow_html=True)
        st.divider()

        # Load users
        if os.path.exists("users.json"):
            with open("users.json", "r") as f:
                USERS = json.load(f)
        else:
            USERS = {}

        menu = ["Login", "Sign Up"]
        choice = st.selectbox("Select Option", menu)

        col1, col2, col3 = st.columns([1,2,1])

        with col2:

            if choice == "Login":

                username = st.text_input("👤 Username")
                password = st.text_input("🔑 Password", type="password")

                if st.button("Login", use_container_width=True):
                    if username in USERS and USERS[username] == password:
                        st.session_state.logged_in = True
                        st.session_state.username = username
                        st.success("✅ Login successful!")
                        st.rerun()
                    else:
                        st.error("❌ Invalid username or password")

            elif choice == "Sign Up":

                new_user = st.text_input("👤 New Username")
                new_password = st.text_input("🔑 New Password", type="password")

                if st.button("Create Account", use_container_width=True):

                    if new_user in USERS:
                        st.warning("⚠️ Username already exists")

                    elif new_user == "" or new_password == "":
                        st.warning("⚠️ Please fill all fields")

                    else:
                        USERS[new_user] = new_password

                        with open("users.json", "w") as f:
                            json.dump(USERS, f)

                        st.success("🎉 Account created successfully!")
                        st.info("Now go to Login and sign in.")

        st.stop()

    # -------- AFTER LOGIN --------

    with st.sidebar:
        st.success(f"Logged in as {st.session_state.username}")
        if st.button("🚪 Logout"):
            st.session_state.logged_in = False
            st.rerun()
    st.title("🧠 Disease Tagger & Insight Tool")
    st.markdown("### Analyze biomedical text and discover therapeutic insights")
    
    # Initialize components
    tagger = DiseaseTagger()
    pubmed = PubMedFetcher()
    
    # Sidebar
    st.sidebar.header("⚙️ Settings")
    max_results = st.sidebar.slider("Max PubMed Results", 1, 10, 3)
    show_fuzzy = st.sidebar.checkbox("Show fuzzy matches", True)
    
    # Main interface
    col1, col2 = st.columns([2, 1])
    
    with col1:
        st.subheader("📝 Input Text")
        user_input = st.text_area(
            "Paste a biomedical sentence or paragraph here:",
            height=150,
            placeholder="Example: This study explores treatment options for diabetes in adults."
        )
        
        # Sample sentences
        if st.button("📋 Use Sample Text"):
            user_input = SAMPLE_SENTENCES[0]
            st.rerun()
    
    with col2:
        st.subheader("📊 Statistics")
        st.metric("Total Diseases", len(DISEASE_AREA_MAP))
        
        # Therapeutic area distribution
        areas = list(set(DISEASE_AREA_MAP.values()))
        area_counts = {area: sum(1 for v in DISEASE_AREA_MAP.values() if v == area) for area in areas}
        st.bar_chart(pd.DataFrame(list(area_counts.items()), columns=['Area', 'Count']).set_index('Area'))
    
    # Process input
    if user_input:
        diseases = tagger.tag_diseases_in_text(user_input)
        
        if diseases:
            st.success(f"✅ Found {len(diseases)} disease(s)")
            
            for i, (disease, area) in enumerate(diseases):
                with st.expander(f"🔬 {disease.title()} ({area})", expanded=True):
                    col1, col2 = st.columns([1, 1])
                    
                    with col1:
                        st.markdown(f"**Disease**: {disease.title()}")
                        st.markdown(f"**Therapeutic Area**: {area}")
                    
                    with col2:
                        if st.button(f"📄 Fetch Abstracts", key=f"fetch_{i}"):
                            with st.spinner("Fetching PubMed insights..."):
                                abstracts = pubmed.fetch_abstracts(disease, max_results)
                                st.text_area(
                                    f"Recent Abstracts for {disease}",
                                    abstracts,
                                    height=300,
                                    key=f"abstracts_{i}"
                                )
        else:
            st.warning("❌ No known diseases detected in the text.")
            
            if show_fuzzy:
                # Try fuzzy matching
                words = user_input.lower().split()
                suggestions = []
                for word in words:
                    similar = tagger.find_similar_diseases(word, n=3, cutoff=0.3)
                    suggestions.extend(similar)
                
                if suggestions:
                    st.info(f"💡 Did you mean: {', '.join(set(suggestions))}?")
    
    # Batch processing
    st.subheader("📦 Batch Processing")
    
    col1, col2 = st.columns([1, 1])
    
    with col1:
        if st.button("🧪 Process Sample Sentences"):
            with st.spinner("Processing..."):
                df = tagger.process_sentences(SAMPLE_SENTENCES)
                st.dataframe(df)
                
                # Download option
                csv = df.to_csv(index=False)
                st.download_button(
                    label="📥 Download Results",
                    data=csv,
                    file_name="tagged_sentences.csv",
                    mime="text/csv"
                )
    
    with col2:
        uploaded_file = st.file_uploader("📁 Upload Text File", type=["txt", "csv"])
        if uploaded_file:
            content = uploaded_file.read().decode("utf-8")
            sentences = content.split('\n')
            sentences = [s.strip() for s in sentences if s.strip()]
            
            if sentences:
                df = tagger.process_sentences(sentences)
                st.dataframe(df)

def main():
    """Main function to run the application."""
    import sys
    
    # Save disease map only once using cache
    save_disease_map_once()
    
    if len(sys.argv) > 1 and sys.argv[1] == "--cli":
        # Run command-line version
        run_command_line_demo()
    else:
        # Run Streamlit app
        create_streamlit_app()

if __name__ == "__main__":
    main()
