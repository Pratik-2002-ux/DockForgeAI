#!/usr/bin/env python3
"""Apply the frozen Track 1 Ki regressor to a prepared supplier library."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.prospective_screening.inference import run_screening


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--model", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    resolve = lambda path: path if path.is_absolute() else ROOT / path
    manifest = run_screening(resolve(args.input), resolve(args.model), resolve(args.output_dir), ROOT, args.overwrite)
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
