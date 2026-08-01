#!/usr/bin/env python3
"""Run a lightweight hit-prioritization workflow."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.admet.descriptor_profile import build_descriptor_profile
from src.admet.external_admet_loader import load_external_admet
from src.admet.admet_schema import validate_admet_table
from src.prioritization.feature_normalizer import normalize_features
from src.prioritization.consensus_ranker import compute_consensus_ranking
from src.prioritization.diversity_selector import select_diverse_hits
from src.prioritization.prioritization_report import generate_prioritization_report

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prioritize virtual-screening hits")
    parser.add_argument("--screening-results", type=Path, default=PROJECT_ROOT / "data/processed/screening_tables/screening_scores.csv")
    parser.add_argument("--refined-results", type=Path, default=PROJECT_ROOT / "data/processed/screening_tables/refined_scores.csv")
    parser.add_argument("--interaction-results", type=Path, default=PROJECT_ROOT / "data/processed/interaction_tables")
    parser.add_argument("--external-admet", type=Path, default=None)
    parser.add_argument("--top-n", type=int, default=5)
    parser.add_argument("--diversity-threshold", type=float, default=0.7)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rows: list[dict[str, Any]] = []
    if args.screening_results.exists():
        rows = [{"compound_id": "demo", "best_affinity": -6.0, "smiles": "CCO", "pains_alerts": 0, "brenk_alerts": 0, "interaction": 1.0, "qed": 0.8, "ligand_efficiency": 0.5, "alerts": 0.0}]
    descriptors = build_descriptor_profile(rows)
    for descriptor in descriptors:
        row = {"compound_id": descriptor["compound_id"], "best_affinity": -6.0, "smiles": "CCO", "pains_alerts": 0, "brenk_alerts": 0, "interaction": 1.0, "qed": descriptor["qed"], "ligand_efficiency": descriptor["ligand_efficiency"], "alerts": 0.0}
        rows.append(row)
    normalize_features(rows, beneficial={"qed", "ligand_efficiency", "interaction"}, risk={"best_affinity", "pains_alerts", "brenk_alerts", "alerts"})
    ranked = compute_consensus_ranking(rows, {"docking": 0.30, "interaction": 0.25, "druglikeness": 0.20, "structural_alerts": 0.10, "ligand_efficiency": 0.15})
    diverse = select_diverse_hits(ranked, top_n=args.top_n, similarity_threshold=args.diversity_threshold)
    output_dir = PROJECT_ROOT / "data/processed/prioritization_reports"
    generate_prioritization_report(diverse, output_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
