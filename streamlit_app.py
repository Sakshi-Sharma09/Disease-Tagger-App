import streamlit as st
import json

with open("disease_map.json", "r") as f:
    disease_area_map = json.load(f)


import streamlit as st
import json
from Bio import Entrez

# ✅ Load disease map once, properly
with open("disease_map.json", "r") as f:
    disease_area_map = json.load(f)

# Email for Entrez API
Entrez.email = "your_email@example.com"

def fetch_pubmed_abstracts(query, max_results=3):
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        record = Entrez.read(handle)
        ids = ",".join(record["IdList"])
        
        fetch = Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="text")
        return fetch.read()
    except Exception as e:
        return f"Error fetching abstracts: {e}"

st.title("🧠 Disease Tagger & Insight Tool")

user_input = st.text_area("Paste a biomedical sentence or paragraph here:")

if user_input:
    matched = False
    for disease in disease_area_map:
        if disease in user_input.lower():
            area = disease_area_map[disease]
            st.markdown(f"### 🔬 Disease Detected: `{disease.title()}`")
            st.markdown(f"**Therapeutic Area**: {area}")
            with st.spinner("Fetching PubMed insights..."):
                abstracts = fetch_pubmed_abstracts(disease)
                st.text_area("📄 Recent Abstracts from PubMed", abstracts, height=300)
            matched = True
            break
    if not matched:
        st.warning("❌ No known disease detected in the text.")
