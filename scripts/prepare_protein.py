#!/usr/bin/env python3
"""Prepare the cached Phase 4 receptor with provenance and reports."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.docking.protein import prepare_receptor


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("data/processed/phase4/protein"))
    parser.add_argument("--pdb-id", default="1HSG")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    output = args.output_dir if args.output_dir.is_absolute() else ROOT / args.output_dir
    print(json.dumps(prepare_receptor(ROOT, output, args.pdb_id, args.overwrite), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
