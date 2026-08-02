#!/usr/bin/env python
"""Build the curated PubChem external-validation dataset."""
from __future__ import annotations
import argparse,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))
from src.external_validation.pubchem_assay_client import PubChemAssayClient
from src.external_validation.pubchem_external_curator import build_external_dataset
def main():
 p=argparse.ArgumentParser(); p.add_argument("--accepted-assays",required=True); p.add_argument("--offline",action="store_true"); p.add_argument("--resume",action="store_true"); p.add_argument("--strict",action="store_true"); p.add_argument("--aid",type=int,nargs="*"); p.add_argument("--refresh",action="store_true"); p.add_argument("--minimum-usable-compounds",type=int,default=30); p.add_argument("--overwrite",action="store_true"); a=p.parse_args(); client=PubChemAssayClient(ROOT/"data/raw/pubchem_external_validation",a.offline); print(build_external_dataset(client,Path(a.accepted_assays),ROOT/"data/processed/pubchem_external_validation/dataset",ROOT/"data/processed/xai_dataset/hivpr_xai_2000_compounds.csv",a.minimum_usable_compounds,a.refresh))
if __name__=="__main__": main()
