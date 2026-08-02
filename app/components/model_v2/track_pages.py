"""Independent Track 1/Track 2 status views."""
import json
from pathlib import Path
import pandas as pd
import streamlit as st
ROOT=Path(__file__).resolve().parents[3]
TRACK_PAGES=["Overview","Track 1 — Ki Dataset","Track 1 — Regression Models","Track 1 — Evaluation","Track 1 — Applicability Domain","Track 1 — Uncertainty","Track 1 — Explainability","Track 1 — Prospective Screening","Track 1 — Candidate Ranking","Track 1 — Reports","Track 2 — Dataset Status","Track 2 — Readiness","Track 2 — Source Balance","Track 2 — Shortcut Analysis","Track 2 — Missing Requirements","Track 2 — Future Binary Model","Track 2 — Reports"]
def _json(path): return json.loads(path.read_text()) if path.exists() else {}
def render(page):
 st.header(f"DockForgeAI v2 — {page}");best=_json(ROOT/"data/model_v2/regression/evaluation/best_model.json");ready=_json(ROOT/"data/model_v2/readiness/readiness_report.json")
 if page=="Overview":st.json({"track_1":"active — Ki-only regression","track_2":"training blocked — data development only"});return
 if page.startswith("Track 2"):
  st.error("TRAINING BLOCKED")
  if page.endswith("Source Balance"):st.json(ready.get("source_counts",{}))
  elif page.endswith("Future Binary Model"):st.code((ROOT/"data/model_v2/binary/README.md").read_text())
  else:st.json({k:ready.get(k) for k in ("training_allowed","final_active_count","final_inactive_count","source_shortcut_risk","blockers","required_next_action")})
  return
 if not best:st.info("Not available: run scripts/model_v2/run_complete_ki_regression.py");return
 if page.endswith("Ki Dataset"):st.dataframe(pd.read_csv(ROOT/"data/model_v2/regression/dataset/ki_regression_dataset.csv"),hide_index=True)
 elif page.endswith("Regression Models"):st.dataframe(pd.read_csv(ROOT/"data/model_v2/regression/evaluation/model_comparison.csv"),hide_index=True)
 elif page.endswith("Evaluation"):st.json(best.get("locked_scaffold_test_metrics"))
 elif page.endswith("Applicability Domain"):st.json(best.get("applicability_domain"))
 elif page.endswith("Uncertainty"):st.json(best.get("uncertainty"))
 elif page.endswith("Explainability"):st.dataframe(pd.read_csv(ROOT/"data/model_v2/regression/explanations/feature_importance.csv").head(30),hide_index=True);st.caption("Model behaviour only; not biological causality.")
 elif page.endswith("Reports"):st.info("Reports are available under reports/model_v2/.")
 else:st.info("Not available: a validated local supplier library is required. Track 2 never participates in prospective screening.")
