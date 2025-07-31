#Disease : Therapeutic Area Map

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

#Sample Biomedical Sentences

sentences = [
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

#Tagging Logic
results = []
for text in sentences:
    tagged = False
    for disease in disease_area_map:
        if disease in text.lower():
            print(f"Text: {text}")
            print(f"Disease: {disease}")
            print(f"Therapeutic Area: {disease_area_map[disease]}\n")
            tagged = True
            break
    if not tagged:
        print(f"Text: {text}")
        print("Disease: Not found\n")

import pandas as pd
df = pd.DataFrame(results)
df

from google.colab import files
uploaded = files.upload()


with open("disease_map.json", "r") as f:
    content = f.read()
    print(repr(content))  # shows special characters and format issues

import spacy
nlp = spacy.load("en_core_web_sm")

text = "The patient was diagnosed with diabetes and hypertension. Prescribed metformin and lisinopril."
doc = nlp(text)

print("Named Entities:")
for ent in doc.ents:
    print(f"{ent.text} ({ent.label_})")


from difflib import get_close_matches

query = "sugar disease"
matches = get_close_matches(query.lower(), disease_area_map.keys(), n=3, cutoff=0.4)

print("Closest Matches:", matches)
if matches:
    print("Therapeutic Area:", disease_area_map[matches[0]])

!pip install Bio
import Bio
from Bio import Entrez
Entrez.email = "your_email@example.com"

def fetch_pubmed_abstracts(query, max_results=3):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    ids = ",".join(record["IdList"])

    fetch = Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="text")
    print(fetch.read())

fetch_pubmed_abstracts("diabetes treatment")

!pip install streamlit

import streamlit as st
import json
from Bio import Entrez

# Load disease map
with open("disease_map.json", "r") as f:
    content = f.read()
    print(repr(content))

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


df.to_csv("tagged_sentences.csv", index=False)
