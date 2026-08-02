#!/usr/bin/env python
from __future__ import annotations
import argparse,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))
from src.models_v2.training import train_vector_model
if __name__=="__main__":
 p=argparse.ArgumentParser(); p.add_argument("--dataset",required=True); p.add_argument("--output",default="data/model_v2/models"); p.add_argument("--model",default="random_forest"); p.add_argument("--seed",type=int,default=2026); a=p.parse_args(); print(train_vector_model(a.dataset,ROOT/a.output,a.model,a.seed))
