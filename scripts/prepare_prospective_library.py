#!/usr/bin/env python3
import argparse,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import pandas as pd
from src.prospective_screening.library_validator import validate_library
if __name__=="__main__":
 p=argparse.ArgumentParser();p.add_argument("--input",required=True);p.add_argument("--output-dir",required=True);a=p.parse_args();d=pd.read_csv(a.input);validate_library(d);o=Path(a.output_dir);o.mkdir(parents=True,exist_ok=True);d.to_csv(o/"candidate_pool.csv",index=False)
