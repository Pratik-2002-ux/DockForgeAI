#!/usr/bin/env python3
"""Validate an existing Phase XAI-1 HIVPR dataset and its companion outputs."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0, str(ROOT))
from src.datasets.dataset_validator import validate_output_files, validate_records  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--dataset", type=Path, required=True); args = parser.parse_args()
    dataset = args.dataset if args.dataset.is_absolute() else ROOT / args.dataset
    records = pd.read_csv(dataset).where(pd.notna(pd.read_csv(dataset)), None).to_dict("records")
    result = validate_records(records)
    stem = dataset.name.removesuffix(".csv"); prefix = stem.removesuffix("_compounds")
    result.update(validate_output_files(dataset, dataset.with_name(f"{prefix}_compounds.parquet"), dataset.with_name(f"{prefix}_compounds.sdf")))
    manifest = json.loads((dataset.parent / "hivpr_xai_dataset_manifest.json").read_text())
    if manifest["final_dataset_size"] != result["row_count"]: raise RuntimeError("Manifest row count mismatch")
    print(json.dumps({"status": "valid", **result}, indent=2)); return 0


if __name__ == "__main__": raise SystemExit(main())
