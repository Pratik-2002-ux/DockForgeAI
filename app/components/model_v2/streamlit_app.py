"""Standalone v2 dashboard; does not import or mutate the v1 app."""
from pathlib import Path
import pandas as pd
import streamlit as st
ROOT=Path(__file__).resolve().parents[3]
st.set_page_config(page_title="DockForgeAI v2",layout="wide")
st.title("DockForgeAI v2")
page=st.sidebar.radio("v2 Navigation",["Model Comparison","Training","Prediction","Explanation","Uncertainty","Calibration","External Validation"])
st.warning("v2 is isolated from the frozen v1 system. Predictions are research outputs and do not establish biological or therapeutic efficacy.")
if page=="Model Comparison":
 path=ROOT/"data/model_v2/model_comparison.csv"; st.dataframe(pd.read_csv(path) if path.exists() else pd.DataFrame())
elif page=="Training": st.info("Training runs only from the explicit v2 CLI; opening this page starts no computation.")
elif page=="Prediction": st.info("Load a versioned v2 artifact after a validated experiment is available.")
elif page=="Explanation": st.info("Explanations describe model behaviour, not causal mechanisms.")
elif page=="Uncertainty": st.info("Ensemble, conformal, and OOD outputs appear here when generated.")
elif page=="Calibration": st.info("ECE, Brier score, temperature scaling, and reliability artifacts appear here.")
else: st.info("Leave-source-out and independent validation artifacts appear here without modifying v1.")
