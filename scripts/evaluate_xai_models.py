#!/usr/bin/env python3
"""Model evaluation is performed without refitting by train_xai_models.py."""
from pathlib import Path
import json
path=Path(__file__).resolve().parents[1]/"data/processed/xai_benchmarks/best_model.json"
if __name__=="__main__": print(json.dumps(json.loads(path.read_text()),indent=2))
