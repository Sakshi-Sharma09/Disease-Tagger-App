import streamlit as st
import json
import re
import pandas as pd
from Bio import Entrez

# Set up the Streamlit app
st.set_page_config(page_title="Disease Tagger App", page_icon="🧬", layout="wide")
st.markdown("<h1 style='text-align: center; color: #4B8BBE;'>🔬 Disease Tagger App</h1>", unsafe_allow_html=True)

# Sidebar with instructions
with st.sidebar:
    st.header("📝 Instructions")
    st.write("""
    1. Paste a biomedical abstract or custom text in the box.
    2. Click 'Tag Diseases' to highlight disease terms.
    3. You can also search PubMed by keyword.
    """)

# Load disease map with caching
@st.cache_data
def load_disease_map():
    with open("disease_map.json", "r") as f:
        return json.load(f)

disease_area_map = load_disease_map()

# Flatten disease terms for quick lookup
disease_terms = []
for terms in disease_area_map.values():
    disease_terms.extend(terms)

# Function to highlight keywords
def highlight_keywords(text, keywords):
    pattern = '|'.join(re.escape(k) for k in keywords if k.strip())
    return re.sub(f"({pattern})", r"<span style='background-color: #feca57;'>\1</span>", text, flags=re.IGNORECASE)

# Fetch PubMed abstract using Biopython
@st.cache_data(ttl=600)
def fetch_pubmed_abstract(query):
    Entrez.email = "your-email@example.com"  # Replace with your email
    handle = Entrez.esearch(db="pubmed", term=query, retmax=1)
    record = Entrez.read(handle)
    handle.close()
    if record["IdList"]:
        article_id = record["IdList"][0]
        fetch_handle = Entrez.efetch(db="pubmed", id=article_id, rettype="abstract", retmode="text")
        abstract_text = fetch_handle.read()
        fetch_handle.close()
        return abstract_text
    else:
        return "No results found."

# Text input
st.subheader("📄 Paste your biomedical text or abstract")
user_input = st.text_area("Enter text here:", height=200)

# PubMed search
st.subheader("🔍 Or fetch a PubMed abstract")
search_query = st.text_input("Enter PubMed keyword:")
if st.button("Fetch from PubMed") and search_query:
    user_input = fetch_pubmed_abstract(search_query)
    st.success("Abstract loaded successfully!")

# Tag and display results
if st.button("🧬 Tag Diseases") and user_input:
    highlighted_text = highlight_keywords(user_input, disease_terms)

    # Extract tagged terms for summary
    tagged_terms = [term for term in disease_terms if re.search(rf"\b{re.escape(term)}\b", user_input, flags=re.IGNORECASE)]

    st.markdown("### 🖍️ Highlighted Text")
    st.markdown(highlighted_text, unsafe_allow_html=True)

    if tagged_terms:
        st.markdown("### 📊 Tagged Disease Terms")
        df = pd.DataFrame({"Disease Term": list(set(tagged_terms))})
        st.dataframe(df, use_container_width=True)

        st.download_button("📥 Download Tagged Text", data=highlighted_text, file_name="tagged_output.html", mime="text/html")
    else:
        st.info("No disease terms found in the input.")
