import streamlit as st
import pandas as pd

# 🔐 Protect Page
if "logged_in" not in st.session_state or not st.session_state.logged_in:
    st.warning("Please login first.")
    st.stop()

st.title("🧪 Pharmacovigilance (PV) Module")

uploaded_file = st.file_uploader("Upload ICSR CSV")

if uploaded_file:
    df = pd.read_csv(uploaded_file)

    st.subheader("Adverse Event Frequency")
    if "Event" in df.columns:
        st.bar_chart(df["Event"].value_counts())

    st.subheader("Seriousness Distribution")
    if "Seriousness" in df.columns:
        st.bar_chart(df["Seriousness"].value_counts())

    st.dataframe(df.head())
