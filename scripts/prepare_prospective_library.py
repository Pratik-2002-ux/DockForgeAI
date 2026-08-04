#!/usr/bin/env python3
"""Prepare a supplier snapshot without inference, docking, or ranking."""
import argparse,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.prospective_screening.library_preparation import prepare_library
if __name__=="__main__":
 p=argparse.ArgumentParser(description=__doc__);p.add_argument("--input",required=True,type=Path);p.add_argument("--output-dir",required=True,type=Path);p.add_argument("--overwrite",action="store_true");p.add_argument("--strict",action="store_true");a=p.parse_args();input_path=a.input if a.input.is_absolute() else ROOT/a.input;output_dir=a.output_dir if a.output_dir.is_absolute() else ROOT/a.output_dir;print(json.dumps(prepare_library(input_path,output_dir,ROOT,a.overwrite,a.strict),indent=2))
