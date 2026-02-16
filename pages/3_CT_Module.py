import streamlit as st
import pandas as pd

# 🔐 Protect Page
if "logged_in" not in st.session_state or not st.session_state.logged_in:
    st.warning("Please login first.")
    st.stop()

st.title("🧬 Clinical Trial (CT) Module")

uploaded_file = st.file_uploader("Upload Clinical Trial CSV")

if uploaded_file:
    df = pd.read_csv(uploaded_file)

    if "Arm" in df.columns and "Event" in df.columns:
        st.subheader("Adverse Events by Trial Arm")
        st.bar_chart(df.groupby("Arm")["Event"].count())

    st.dataframe(df.head())
