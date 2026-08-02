"""Transparent multi-criterion model ranking."""
from __future__ import annotations
import pandas as pd
DEFAULT_WEIGHTS={"roc_auc":.15,"pr_auc":.2,"mcc":.2,"balanced_accuracy":.15,"calibration":.1,"external_validation":.1,"novel_scaffold":.05,"inference_speed":.025,"memory_usage":.025}
def rank_models(frame,weights=None):
 weights=weights or DEFAULT_WEIGHTS; data=frame.copy(); score=0.
 for column,weight in weights.items():
  values=data[column].astype(float); normalized=(values-values.min())/(values.max()-values.min()) if values.max()>values.min() else values*0+1; normalized=1-normalized if column in {"calibration","inference_speed","memory_usage"} else normalized; score=score+weight*normalized
 data["selection_score"]=score; return data.sort_values("selection_score",ascending=False).reset_index(drop=True)
