"""Read-only prospective-discovery status pages."""
import json
from pathlib import Path
import streamlit as st
ROOT=Path(__file__).resolve().parents[3]
PAGES=["Discovery overview","Ki regression","Prospective library","ML screening","Applicability domain","Candidate diversity","Docking and interactions","Frozen shortlist","Biochemical study","Biochemical results","Orthogonal evidence","Cellular evidence","ADME and safety","Evidence matrix","Candidate status","Reports","Reproducibility"]
def render(page):
 st.header(f"DockForgeAI Prospective Discovery — {page}");st.warning("Computational predictions are not experimental evidence and do not establish inhibition or therapeutic activity.")
 best=ROOT/"data/model_v2/regression/evaluation/best_model.json";shortlists=list((ROOT/"data/processed/prospective_studies").glob("*/shortlist_manifest.json")) if (ROOT/"data/processed/prospective_studies").exists() else []
 if page=="Ki regression" and best.exists():st.json(json.loads(best.read_text()));return
 if page in {"Discovery overview","Reproducibility"}:st.json({"phase_1_ki_regression":best.exists(),"phase_2_frozen_shortlist":bool(shortlists),"phase_3_biochemical_evidence":False,"phase_4_orthogonal_evidence":False,"phase_5_cellular_adme_safety":False});return
 requirements={"Prospective library":"A local CSV with compound_id, canonical_smiles, and supplier_id.","ML screening":"A validated local library and frozen Phase-1 model.","Frozen shortlist":"A completed computational review followed by an immutable freeze command.","Biochemical study":"A frozen shortlist and qualified-laboratory preregistration.","Biochemical results":"Genuine collaborator results matching the registered schema.","Orthogonal evidence":"QC-passing independent collaborator evidence.","Cellular evidence":"QC-passing qualified-laboratory cellular results.","ADME and safety":"QC-passing collaborator ADME/safety results."}
 st.info("Not available");st.code(requirements.get(page,"Required phase artifacts have not been supplied. No result is fabricated."))
