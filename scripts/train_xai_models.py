#!/usr/bin/env python3
from __future__ import annotations
import argparse,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.ml.experiment_runner import FEATURES,MODELS,run_experiments
def main():
 p=argparse.ArgumentParser();p.add_argument("--all-models",action="store_true");p.add_argument("--all-feature-sets",action="store_true");p.add_argument("--seeds",nargs="+",type=int,default=[2026]);p.add_argument("--overwrite",action="store_true");a=p.parse_args();print(run_experiments(ROOT/"data/processed/xai_features",ROOT/"data/processed/xai_models",ROOT/"data/processed/xai_predictions",ROOT/"data/processed/xai_benchmarks",a.seeds));return 0
if __name__=="__main__":raise SystemExit(main())
