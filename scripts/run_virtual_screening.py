#!/usr/bin/env python3
"""Run a lightweight virtual-screening demonstration workflow."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.ligands.library_loader import load_library
from src.ligands.standardizer import standardize_ligands
from src.ligands.druglikeness import calculate_descriptors
from src.ligands.conformer_generator import generate_conformers
from src.ligands.pdbqt_batch import prepare_pdbqt
from src.screening.screening_config import write_configs
from src.screening.batch_runner import run_batch
from src.screening.score_parser import parse_scores
from src.screening.hit_selector import select_hits
from src.screening.refinement_runner import run_refinement
from src.screening.screening_report import generate_screening_report

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run a lightweight virtual-screening demonstration")
    parser.add_argument("--input", type=Path, default=PROJECT_ROOT / "data/raw/screening_libraries/demo_library.csv")
    parser.add_argument("--max-compounds", type=int, default=10)
    parser.add_argument("--demo-library", action="store_true")
    parser.add_argument("--refine", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.demo_library:
        demo_dir = PROJECT_ROOT / "data/raw/screening_libraries"
        demo_dir.mkdir(parents=True, exist_ok=True)
        demo_file = demo_dir / "demo_library.csv"
        if not demo_file.exists():
            demo_file.write_text("compound_id,compound_name,smiles,source\n1,DemoA,CCO,demo\n2,DemoB,CCN,demo\n3,DemoC,CC(C)O,demo\n", encoding="utf-8")
        args.input = demo_file

    library = load_library(args.input, max_compounds=args.max_compounds)
    standardized = standardize_ligands(library["records"])
    descriptors = calculate_descriptors(standardized["standardized"])
    conformers = generate_conformers(standardized["standardized"], PROJECT_ROOT / "data/processed/screening_library/sdf")
    pdbqt = prepare_pdbqt(conformers, PROJECT_ROOT / "data/processed/screening_library/pdbqt", overwrite=args.overwrite)
    configs = write_configs(PROJECT_ROOT / "data/processed/screening_runs")
    batch_results = run_batch(standardized["standardized"], configs["screening_config"], PROJECT_ROOT / "data/processed/screening_runs", timeout_seconds=1800, resume=args.resume)
    score_rows = parse_scores(batch_results["results"], PROJECT_ROOT / "data/processed/screening_tables/screening_scores.csv")
    if args.refine:
        selected_hits = select_hits(PROJECT_ROOT / "data/processed/screening_tables/screening_scores.csv", PROJECT_ROOT / "data/processed/screening_tables/refinement_candidates.csv")
        refined_rows = run_refinement(selected_hits, PROJECT_ROOT / "data/processed/screening_tables/refined_scores.csv")
    else:
        refined_rows = []
    generate_screening_report({"library_size": len(standardized["standardized"]), "accepted_count": len(standardized["standardized"]), "rejected_count": len(standardized["rejected"]), "screening_result_count": len(score_rows), "refined_count": len(refined_rows)}, PROJECT_ROOT / "data/processed/screening_reports")
    return 0


if __name__ == "__main__":
    sys.exit(main())
