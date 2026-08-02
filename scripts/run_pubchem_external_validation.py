#!/usr/bin/env python
"""Evaluate the locked model on the curated PubChem dataset."""
from __future__ import annotations
import argparse,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))
from src.external_validation.pubchem_external_evaluator import evaluate_locked_model
from src.external_validation.pubchem_external_report import generate_report
def main():
 p=argparse.ArgumentParser(); p.add_argument("--locked-model",required=True); p.add_argument("--bootstrap-iterations",type=int,default=1000); p.add_argument("--offline",action="store_true"); p.add_argument("--resume",action="store_true"); p.add_argument("--strict",action="store_true"); p.add_argument("--aid",type=int,nargs="*"); p.add_argument("--refresh",action="store_true"); p.add_argument("--overwrite",action="store_true"); a=p.parse_args(); root=ROOT/"data/processed/pubchem_external_validation"; print(evaluate_locked_model(root/"dataset/pubchem_external_validation.csv",root/"evaluation",a.bootstrap_iterations,a.overwrite)); generate_report(root)
if __name__=="__main__": main()
