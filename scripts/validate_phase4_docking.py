#!/usr/bin/env python3
"""Generate Phase 4 reports from computed crystallographic redocking validation."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.docking.validation import validate_docking
def main()->int:
 p=argparse.ArgumentParser();p.add_argument("--output-dir",type=Path,default=Path("data/processed/phase4/validation"));p.add_argument("--overwrite",action="store_true");a=p.parse_args();o=a.output_dir if a.output_dir.is_absolute() else ROOT/a.output_dir;print(json.dumps(validate_docking(ROOT/"data/processed/docking_results/1HSG_MK1_redocking_summary.json",ROOT/"data/processed/docking_tables/1HSG_MK1_redocking_analysis.csv",o,ROOT,a.overwrite),indent=2));return 0
if __name__=="__main__":raise SystemExit(main())
