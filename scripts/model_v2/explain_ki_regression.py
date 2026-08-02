#!/usr/bin/env python3
"""Generate conservative feature importance for the frozen Ki model."""
import json
from pathlib import Path
import joblib,numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[2]
if __name__=="__main__":
 m=joblib.load(ROOT/"data/model_v2/regression/models/best_ki_regressor.joblib");out=ROOT/"data/model_v2/regression/explanations";out.mkdir(parents=True,exist_ok=True);values=getattr(m,"feature_importances_",np.array([]));pd.DataFrame({"feature_index":range(len(values)),"importance":values}).sort_values("importance",ascending=False).to_csv(out/"feature_importance.csv",index=False);(out/"explanation_manifest.json").write_text(json.dumps({"describes_model_behavior_only":True,"biological_causality":False,"hash_collisions_possible":True},indent=2)+"\n")
