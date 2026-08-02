#!/usr/bin/env python3
import argparse,json,sys,pandas as pd
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.prospective_screening.shortlist_freezer import freeze_shortlist
if __name__=="__main__":
 p=argparse.ArgumentParser();p.add_argument("--input",required=True);p.add_argument("--output-dir",required=True);p.add_argument("--metadata",required=True);a=p.parse_args();freeze_shortlist(pd.read_csv(a.input),a.output_dir,json.load(open(a.metadata)))
