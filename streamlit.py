# streamlit_app.py

import streamlit as st
import re
import requests
from typing import List, Dict, Set

# Page configuration
st.set_page_config(
    page_title="Disease Tagger App",
    page_icon="🔬",
    layout="wide"
)

# Title with icon
st.title("🔬 Disease Tagger App")

# Instructions sidebar
with st.sidebar:
    st.header("📝 Instructions")
    st.write("1. Paste a biomedical abstract or custom text in the box.")
    st.write("2. Click 'Tag Diseases' to highlight disease terms.")
    st.write("3. You can also search PubMed by keyword.")

# Disease terms database (expanded list)
DISEASE_TERMS = {
    # Cardiovascular
    "hypertension", "stroke", "heart disease", "myocardial infarction", "arrhythmia", 
    "congestive heart failure", "atherosclerosis", "coronary artery disease",
    
    # Neurological
    "alzheimer's disease", "parkinson's disease", "epilepsy", "multiple sclerosis", 
    "migraine", "dementia", "stroke", "seizure", "neuropathy",
    
    # Cancer/Oncology
    "cancer", "tumor", "carcinoma", "sarcoma", "leukemia", "lymphoma", 
    "breast cancer", "lung cancer", "prostate cancer", "melanoma", "glioblastoma",
    
    # Endocrine
    "diabetes", "hypothyroidism", "hyperthyroidism", "cushing's syndrome", 
    "acromegaly", "metabolic syndrome", "obesity",
    
    # Respiratory
    "asthma", "copd", "chronic obstructive pulmonary disease", "tuberculosis", 
    "pneumonia", "pulmonary fibrosis", "lung disease",
    
    # Infectious
    "covid-19", "malaria", "dengue", "hiv", "aids", "hepatitis", "influenza",
    "tuberculosis", "sepsis", "infection",
    
    # Autoimmune/Inflammatory
    "rheumatoid arthritis", "lupus", "crohn's disease", "ulcerative colitis", 
    "psoriasis", "multiple sclerosis", "inflammatory bowel disease",
    
    # Kidney
    "chronic kidney disease", "renal failure", "glomerulonephritis", "nephritis",
    
    # Mental Health
    "depression", "anxiety", "schizophrenia", "bipolar disorder", "ptsd",
    
    # Other common terms
    "ischemia", "thrombosis", "embolism", "fibrosis", "cirrhosis", "osteoporosis",
    "glaucoma", "macular degeneration", "cataracts"
}

# Therapeutic area mapping
DISEASE_AREA_MAP = {
    "hypertension": "Cardiology",
    "stroke": "Cardiology/Neurology", 
    "heart disease": "Cardiology",
    "chronic kidney disease": "Nephrology",
    "diabetes": "Endocrinology",
    "asthma": "Pulmonology",
    "cancer": "Oncology",
    "alzheimer's disease": "Neurology",
    "depression": "Psychiatry",
    # Add more mappings as needed
}

def extract_diseases_from_text(text: str) -> List[Dict]:
    """Extract disease terms from text and return with positions"""
    found_diseases = []
    text_lower = text.lower()
    
    # Sort disease terms by length (longest first) to avoid partial matches
    sorted_diseases = sorted(DISEASE_TERMS, key=len, reverse=True)
    
    for disease in sorted_diseases:
        # Use word boundaries to avoid partial matches
        pattern = r'\b' + re.escape(disease.lower()) + r'\b'
        matches = re.finditer(pattern, text_lower)
        
        for match in matches:
            # Get the original case from the text
            original_text = text[match.start():match.end()]
            found_diseases.append({
                'term': original_text,
                'start': match.start(),
                'end': match.end(),
                'therapeutic_area': DISEASE_AREA_MAP.get(disease.lower(), "Unknown")
            })
    
    # Remove duplicates and overlapping matches
    found_diseases = remove_overlapping_matches(found_diseases)
    return found_diseases

def remove_overlapping_matches(matches: List[Dict]) -> List[Dict]:
    """Remove overlapping matches, keeping the longest ones"""
    if not matches:
        return matches
    
    # Sort by start position
    matches.sort(key=lambda x: x['start'])
    
    filtered = []
    for match in matches:
        # Check if this match overlaps with any already added
        overlaps = False
        for existing in filtered:
            if (match['start'] < existing['end'] and match['end'] > existing['start']):
                overlaps = True
                break
        
        if not overlaps:
            filtered.append(match)
    
    return filtered

def highlight_diseases_in_text(text: str, diseases: List[Dict]) -> str:
    """Highlight disease terms in the text using HTML"""
    if not diseases:
        return text
    
    # Sort by start position in reverse order to maintain positions when inserting HTML
    diseases_sorted = sorted(diseases, key=lambda x: x['start'], reverse=True)
    
    highlighted_text = text
    for disease in diseases_sorted:
        start, end = disease['start'], disease['end']
        term = disease['term']
        therapeutic_area = disease['therapeutic_area']
        
        # Create highlighted span with tooltip
        highlighted_term = f'<span style="background-color: #ffeb3b; padding: 2px 4px; border-radius: 3px; font-weight: bold;" title="Therapeutic Area: {therapeutic_area}">{term}</span>'
        
        highlighted_text = highlighted_text[:start] + highlighted_term + highlighted_text[end:]
    
    return highlighted_text

def search_pubmed_mock(keyword: str) -> str:
    """Mock PubMed search - returns a sample abstract"""
    sample_abstracts = {
        "hypertension": "Hypertension is a major risk factor for stroke, heart disease, and chronic kidney disease. Recent studies have shown that early intervention with antihypertensive medications can significantly reduce the risk of cardiovascular events. The prevalence of hypertension continues to increase globally, particularly in developing countries.",
        
        "diabetes": "Type 2 diabetes mellitus affects millions worldwide and is associated with serious complications including diabetic nephropathy, retinopathy, and cardiovascular disease. Management typically involves lifestyle modifications and pharmacological interventions to maintain glycemic control.",
        
        "cancer": "Cancer remains one of the leading causes of mortality worldwide. Recent advances in immunotherapy have shown promising results in treating various types of cancer including melanoma, lung cancer, and breast cancer. Early detection and personalized treatment approaches are crucial for improving patient outcomes."
    }
    
    return sample_abstracts.get(keyword.lower(), f"Sample abstract for {keyword}: This is a mock abstract that would contain information about {keyword} and related medical conditions. In a real implementation, this would fetch actual abstracts from PubMed.")

# Main interface
st.markdown("### 📄 Paste your biomedical text or abstract")

# Text input area
text_input = st.text_area(
    "Enter text here:",
    value="Hypertension is a major risk factor for stroke, heart disease, and chronic kidney disease.",
    height=200,
    placeholder="Paste your biomedical abstract or custom text here..."
)

# PubMed search section
st.markdown("### 🔍 Or fetch a PubMed abstract")
col1, col2 = st.columns([3, 1])

with col1:
    pubmed_keyword = st.text_input("Enter PubMed keyword:", placeholder="e.g., hypertension, diabetes, cancer")

with col2:
    st.write("")  # Empty space for alignment
    fetch_pubmed = st.button("Fetch from PubMed", type="secondary")

# Handle PubMed fetch
if fetch_pubmed and pubmed_keyword:
    mock_abstract = search_pubmed_mock(pubmed_keyword)
    text_input = mock_abstract
    st.success(f"✅ Fetched mock abstract for: {pubmed_keyword}")

# Tag diseases button
if st.button("🏷️ Tag Diseases", type="primary"):
    if text_input.strip():
        # Extract diseases
        found_diseases = extract_diseases_from_text(text_input)
        
        if found_diseases:
            st.markdown("### 🎯 Disease Terms Found:")
            
            # Display highlighted text
            highlighted_text = highlight_diseases_in_text(text_input, found_diseases)
            st.markdown(highlighted_text, unsafe_allow_html=True)
            
            # Display disease summary
            st.markdown("### 📊 Disease Summary:")
            
            # Create a summary table
            disease_summary = []
            therapeutic_areas = set()
            
            for disease in found_diseases:
                disease_summary.append({
                    "Disease Term": disease['term'],
                    "Therapeutic Area": disease['therapeutic_area']
                })
                therapeutic_areas.add(disease['therapeutic_area'])
            
            # Display as a table
            if disease_summary:
                import pandas as pd
                df = pd.DataFrame(disease_summary)
                st.dataframe(df, use_container_width=True)
                
                # Show therapeutic areas count
                st.markdown(f"**Found {len(found_diseases)} disease terms across {len(therapeutic_areas)} therapeutic areas**")
                
                # Show therapeutic areas
                if therapeutic_areas:
                    areas_list = ", ".join(sorted(therapeutic_areas))
                    st.info(f"🏥 Therapeutic Areas: {areas_list}")
        else:
            st.warning("⚠️ No disease terms found in the provided text.")
    else:
        st.error("❌ Please enter some text to analyze.")

# Footer
st.markdown("---")
st.markdown("💡 **Tip:** This tool identifies common disease terms in biomedical text. The therapeutic area mapping helps categorize diseases by medical specialty.")