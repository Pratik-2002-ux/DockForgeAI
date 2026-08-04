#!/usr/bin/env python3
"""Prepare the committed Phase 4 Top-30 ligand review pool."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.docking.ligands import prepare_ligands


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("data/processed/prospective_studies/HIVPR_2026_001/ml_screening/top_30_diverse_review_pool.csv"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/processed/phase4/ligands"))
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(); resolve = lambda path: path if path.is_absolute() else ROOT / path
    print(json.dumps(prepare_ligands(resolve(args.input), resolve(args.output_dir), ROOT, args.overwrite, args.seed), indent=2))
    return 0


if __name__ == "__main__": raise SystemExit(main())
