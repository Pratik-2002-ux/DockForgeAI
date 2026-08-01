#!/usr/bin/env python3
from __future__ import annotations
import argparse,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
from src.features.feature_matrix import build_feature_artifacts
def main()->int:
 p=argparse.ArgumentParser(); p.add_argument("--dataset",type=Path,required=True); p.add_argument("--overwrite",action="store_true"); a=p.parse_args(); dataset=a.dataset if a.dataset.is_absolute() else ROOT/a.dataset
 print(build_feature_artifacts(dataset,ROOT/"data/processed/xai_features",a.overwrite)); return 0
if __name__=="__main__": raise SystemExit(main())
