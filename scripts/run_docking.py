#!/usr/bin/env python3
"""Run real AutoDock Vina docking for the Phase 4 ligand set."""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from src.docking.phase4_runner import run_docking


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ligand-dir", type=Path, default=Path("data/processed/phase4/ligands/pdbqt")); parser.add_argument("--receptor", type=Path, default=Path("data/processed/phase4/protein/receptor_prepared.pdbqt")); parser.add_argument("--box", type=Path, default=Path("data/processed/docking_boxes/1HSG_MK1_box.json")); parser.add_argument("--output-dir", type=Path, default=Path("data/processed/phase4/docking")); parser.add_argument("--engine", default="vina", choices=["vina", "autodock-gpu", "quickvina2"]); parser.add_argument("--seed", type=int, default=2026); parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(); resolve = lambda path: path if path.is_absolute() else ROOT / path
    print(json.dumps(run_docking(resolve(args.ligand_dir), resolve(args.receptor), resolve(args.box), resolve(args.output_dir), ROOT, args.overwrite, args.engine, args.seed), indent=2)); return 0


if __name__ == "__main__": raise SystemExit(main())
