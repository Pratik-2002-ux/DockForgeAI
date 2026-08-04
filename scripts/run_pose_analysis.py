#!/usr/bin/env python3
"""Analyze computed Phase 4 best poses."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.docking.analysis import analyze_poses
def main()->int:
 p=argparse.ArgumentParser();p.add_argument("--receptor",type=Path,default=Path("data/processed/phase4/protein/receptor_prepared.pdb"));p.add_argument("--docking-dir",type=Path,default=Path("data/processed/phase4/docking"));p.add_argument("--output-dir",type=Path,default=Path("data/processed/phase4/analysis"));p.add_argument("--overwrite",action="store_true");a=p.parse_args();r=lambda x:x if x.is_absolute() else ROOT/x;print(json.dumps(analyze_poses(r(a.receptor),r(a.docking_dir),r(a.output_dir),ROOT,a.overwrite),indent=2));return 0
if __name__=="__main__":raise SystemExit(main())
