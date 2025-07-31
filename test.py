# streamlit_app.py

import streamlit as st

# Title
st.title("🧠 Disease to Therapeutic Area Mapper")

# Subtitle
st.write("Enter a disease name below to find out its therapeutic area.")

# Disease-to-Therapeutic-Area Mapping
disease_area_map = {
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

    # Infectious Disease
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


def get_therapeutic_area(disease):
    disease = disease.strip().lower()
    return disease_area_map.get(disease, "❌ Unknown or not mapped")


# Input Box
disease_input = st.text_input("🔍 Enter disease name:")

# Display result
if disease_input:
    area = get_therapeutic_area(disease_input)
    st.success(f"🏥 Therapeutic Area: **{area}**")
