#!/usr/bin/env python
"""Discover and audit PubChem BioAssays using official APIs."""
from __future__ import annotations
import argparse,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))
from src.external_validation.pubchem_assay_client import PubChemAssayClient
from src.external_validation.pubchem_assay_discovery import discover_assays
def main():
 p=argparse.ArgumentParser(); p.add_argument("--query",default="HIV-1 protease"); p.add_argument("--exclude-depositor",default="ChEMBL"); p.add_argument("--offline",action="store_true"); p.add_argument("--resume",action="store_true"); p.add_argument("--strict",action="store_true"); p.add_argument("--aid",type=int,nargs="*"); p.add_argument("--refresh",action="store_true"); p.add_argument("--max-assays",type=int,default=250); p.add_argument("--minimum-usable-compounds",type=int,default=30); p.add_argument("--overwrite",action="store_true"); a=p.parse_args(); out=ROOT/"data/processed/pubchem_external_validation/assay_discovery"; client=PubChemAssayClient(ROOT/"data/raw/pubchem_external_validation",a.offline); manifest,_=discover_assays(client,out,a.query,a.exclude_depositor,a.max_assays,a.minimum_usable_compounds,a.refresh); print(manifest)
if __name__=="__main__": main()
