import streamlit as st
import pandas as pd

# 🔐 Protect Page
if "logged_in" not in st.session_state or not st.session_state.logged_in:
    st.warning("Please login first.")
    st.stop()

st.title("💊 BA/BE Analytics Module")

uploaded_file = st.file_uploader("Upload PK Dataset CSV")

if uploaded_file:
    df = pd.read_csv(uploaded_file)

    if "Cmax" in df.columns and "AUC" in df.columns:
        st.subheader("Cmax & AUC Trends")
        st.line_chart(df[["Cmax", "AUC"]])

    st.dataframe(df.describe())
