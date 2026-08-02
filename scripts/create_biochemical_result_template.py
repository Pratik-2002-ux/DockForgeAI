#!/usr/bin/env python3
import argparse,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import pandas as pd
from src.experimental_validation.biochemical_result_schema import REQUIRED
if __name__=="__main__":
 p=argparse.ArgumentParser();p.add_argument("--output",required=True);a=p.parse_args();Path(a.output).parent.mkdir(parents=True,exist_ok=True);pd.DataFrame(columns=REQUIRED).to_csv(a.output,index=False)
