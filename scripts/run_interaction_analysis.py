#!/usr/bin/env python3
"""Run the Phase 3 protein–ligand interaction analysis workflow."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.analysis.complex_builder import build_complex
from src.analysis.interaction_engine import analyze_complex
from src.analysis.interaction_parser import parse_interactions
from src.analysis.interaction_report import generate_interaction_report

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run interaction analysis")
    parser.add_argument("--pose", choices=["crystal", "rank1", "best-rmsd", "all"], default="all")
    parser.add_argument("--engine", choices=["auto", "plip", "fallback"], default="auto")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    receptor_path = PROJECT_ROOT / "data/processed/protonated_receptors/1HSG_chains-A-B_receptor.pdb"
    ligand_path = PROJECT_ROOT / "data/processed/crystallographic_ligands/1HSG_MK1_crystallographic_corrected.sdf"
    output_dir = PROJECT_ROOT / "data/processed/interaction_results"
    output_dir.mkdir(parents=True, exist_ok=True)

    poses = ["crystal", "rank1", "best-rmsd"] if args.pose == "all" else [args.pose]
    rows: list[dict[str, Any]] = []
    for pose_label in poses:
        complex_path = output_dir / f"{pose_label}_complex.pdb"
        build_complex(receptor_path=receptor_path, ligand_path=ligand_path, output_path=complex_path, label=pose_label)
        result = analyze_complex(complex_path=complex_path, pose_label=pose_label, pose_rank=1 if pose_label == "crystal" else 1, engine=args.engine)
        parsed = parse_interactions(result["interactions"], output_dir=output_dir, pose_label=pose_label, pose_rank=1)
        rows.extend(result["interactions"])
        (output_dir / f"{pose_label}_analysis.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    generate_interaction_report(rows=rows, output_dir=PROJECT_ROOT / "data/processed/interaction_reports")
    return 0


if __name__ == "__main__":
    sys.exit(main())
