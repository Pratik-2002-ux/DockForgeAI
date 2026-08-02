#!/usr/bin/env python
"""Build a v2 dataset from user-supplied experimental source exports."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))
from src.models_v2.curation import curate_records
from src.models_v2.splits import scaffold_split
import pandas as pd
def main():
 p=argparse.ArgumentParser(); p.add_argument("--input",required=True); p.add_argument("--output",default="data/model_v2/dataset.csv"); p.add_argument("--seed",type=int,default=2026); a=p.parse_args(); raw=pd.read_csv(a.input); data,excluded=curate_records(raw); data["dataset_split"]=scaffold_split(data,a.seed); out=ROOT/a.output; out.parent.mkdir(parents=True,exist_ok=True); data.to_csv(out,index=False); excluded.to_csv(out.with_name("exclusions.csv"),index=False); print(json.dumps({"rows":len(data),"excluded":len(excluded)}))
if __name__=="__main__": main()
