"""Read-only v2 status pages for the combined DockForgeAI application."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parents[3]
READINESS = ROOT / "data/model_v2/readiness/readiness_report.json"
DATASET = ROOT / "data/model_v2/dataset/v2_experimental_dataset.csv"


def _readiness() -> dict:
    if not READINESS.is_file():
        return {}
    return json.loads(READINESS.read_text(encoding="utf-8"))


def render_v2_page(page: str) -> None:
    """Render a v2 page without starting training or network operations."""
    st.header(f"DockForgeAI v2 — {page}")
    report = _readiness()
    st.warning(
        "v2 is experimental. Predictions do not establish binding, therapeutic "
        "activity, or biological causality. The locked external benchmark is "
        "never used for training or model selection."
    )
    if not report:
        st.error(f"Missing artifact: {READINESS.relative_to(ROOT)}")
        st.code("conda run -n dockforge-v2 python scripts/model_v2/audit_v2_readiness.py")
        return
    if page in {"Overview", "System Readiness"}:
        st.json({
            "status": report["status"],
            "training_allowed": report["training_allowed"],
            "approved_training_sources": report["approved_training_sources"],
        })
        for blocker in report.get("blockers", []):
            st.error(blocker)
        st.caption(report.get("required_next_input", ""))
        return
    if page == "Dataset":
        st.json(report.get("chembl_preliminary_audit", {}))
        if DATASET.is_file():
            st.dataframe(pd.read_csv(DATASET), hide_index=True)
        else:
            st.info("No v2 training dataset was generated because readiness failed.")
        return
    artifact_map = {
        "Feature Engineering": "data/model_v2/features/feature_manifest.json",
        "Model Training": "data/model_v2/models/model_registry.json",
        "Model Comparison": "data/model_v2/benchmarks/model_comparison.csv",
        "Explainability": "data/model_v2/explanations/explanation_manifest.json",
        "Uncertainty": "data/model_v2/benchmarks/uncertainty_comparison.csv",
        "Calibration": "data/model_v2/benchmarks/calibration_comparison.csv",
        "Applicability Domain": "data/model_v2/benchmarks/applicability_domain.csv",
        "Locked External Validation": "data/model_v2/benchmarks/external_comparison.csv",
        "v1 vs v2 Comparison": "data/model_v2/benchmarks/external_comparison.csv",
        "Reports": "reports/model_v2/v2_report_manifest.json",
    }
    artifact = ROOT / artifact_map[page]
    if artifact.exists():
        st.success(f"Artifact available: {artifact.relative_to(ROOT)}")
    else:
        st.info(f"Not generated: {artifact.relative_to(ROOT)}")
        st.code("Resolve readiness blockers, then run scripts/model_v2/run_complete_v2_pipeline.py --offline --resume")

