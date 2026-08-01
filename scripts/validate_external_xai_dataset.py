#!/usr/bin/env python
"""Run locked-model inference against a user-supplied dataset."""
from __future__ import annotations
import argparse,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))
from src.validation.external_dataset_validator import validate_external_dataset
def main():
    p=argparse.ArgumentParser(); p.add_argument("--input",required=True); p.add_argument("--id-column",required=True); p.add_argument("--smiles-column",required=True); p.add_argument("--label-column"); p.add_argument("--output-dir",required=True); p.add_argument("--strict",action="store_true"); p.add_argument("--allow-scaffold-overlap",action="store_true"); p.add_argument("--overwrite",action="store_true"); a=p.parse_args(); validate_external_dataset(a.input,a.id_column,a.smiles_column,a.output_dir,a.label_column,a.strict,a.allow_scaffold_overlap,a.overwrite)
if __name__=="__main__": main()
