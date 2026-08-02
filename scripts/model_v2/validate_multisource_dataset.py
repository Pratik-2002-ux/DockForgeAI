#!/usr/bin/env python3
import argparse,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from src.models_v2.data.dataset_validator import validate_dataset
if __name__=="__main__":
 p=argparse.ArgumentParser(); p.add_argument("--dataset",required=True); a=p.parse_args(); print(json.dumps(validate_dataset(ROOT/a.dataset),indent=2))
