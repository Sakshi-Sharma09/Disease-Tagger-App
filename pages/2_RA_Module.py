import streamlit as st
import pandas as pd

# 🔐 Protect Page
if "logged_in" not in st.session_state or not st.session_state.logged_in:
    st.warning("Please login first.")
    st.stop()

st.title("📜 Regulatory Affairs (RA) Module")

uploaded_file = st.file_uploader("Upload Submission Tracker CSV")

if uploaded_file:
    df = pd.read_csv(uploaded_file)

    if "Status" in df.columns:
        st.subheader("Submission Status Overview")
        st.bar_chart(df["Status"].value_counts())

    st.dataframe(df)
