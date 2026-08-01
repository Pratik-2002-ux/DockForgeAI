"""Reuse existing repeated-seed predictions for stability analysis."""
from __future__ import annotations
import json, numpy as np, pandas as pd
from scipy.stats import spearmanr
from .core import OUTPUT, PREDICTIONS, BENCHMARKS, BEST_MODEL_ID, write_json

def run_stability(output=OUTPUT/"stability"):
    comparison=pd.read_csv(BENCHMARKS/"model_comparison.csv"); comparison.to_csv(output/"seed_metric_stability.csv",index=False)
    pred=pd.read_csv(PREDICTIONS/"test_predictions.csv"); family=pred[(pred.model_id.str.startswith("logistic_regression__fingerprint"))].copy(); wide=family.pivot(index="compound_id",columns="seed",values="predicted_probability"); summary=wide.agg(["mean","std","min","max"],axis=1).reset_index(); summary["class_assignment_consistency"]=(wide.ge(.5).sum(axis=1).clip(lower=0).map(lambda n:max(n,wide.shape[1]-n))/wide.shape[1]).to_numpy(); summary.to_csv(output/"prediction_stability.csv",index=False); summary[["compound_id","class_assignment_consistency"]].to_csv(output/"classification_consistency.csv",index=False)
    corr=wide.corr(method="spearman"); corr.to_csv(output/"feature_importance_stability.csv")
    result={"seed_count":int(wide.shape[1]),"mean_probability_std":float(summary["std"].mean()),"mean_rank_correlation":float(corr.where(~np.eye(len(corr),dtype=bool)).stack().mean()),"fully_consistent_fraction":float(summary.class_assignment_consistency.eq(1).mean())}; write_json(output/"stability_summary.json",result); return result
