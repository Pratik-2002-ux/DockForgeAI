#!/usr/bin/env python3
import argparse,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.experimental_validation.study_registration import register
if __name__=="__main__":
 p=argparse.ArgumentParser();p.add_argument("--registration-json",required=True);p.add_argument("--output",required=True);a=p.parse_args();register(json.load(open(a.registration_json)),a.output)
