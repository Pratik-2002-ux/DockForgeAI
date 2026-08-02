#!/usr/bin/env python3
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
if __name__=="__main__":print((ROOT/"data/model_v2/regression/evaluation/best_model.json").read_text())
