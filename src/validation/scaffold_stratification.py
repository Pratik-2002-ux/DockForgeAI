"""Scaffold-frequency robustness analysis."""
from __future__ import annotations
import numpy as np, pandas as pd
from .core import OUTPUT, best_predictions, load_data, safe_metrics, write_json

def scaffold_category(n): return "singleton" if n==1 else "low" if n<=3 else "medium" if n<=10 else "high"
def run_scaffold_stratification(output=OUTPUT/"scaffolds"):
    d=load_data(); test=d[d.dataset_split.eq("test")].copy(); p=best_predictions().set_index("compound_id"); test["probability"]=test.compound_id.map(p.predicted_probability); test["category"]=test.scaffold_frequency.map(scaffold_category); train_scaff=set(d.loc[d.dataset_split.eq("train"),"scaffold_id"]); test["seen_in_train"]=test.scaffold_id.isin(train_scaff); class_n=d.groupby("scaffold_id").binary_label.nunique(); test["label_structure"]=test.scaffold_id.map(class_n).map(lambda n:"mixed_labels" if n>1 else "one_class_only")
    rows=[]
    for dimension in ("category","seen_in_train","label_structure"):
        for name,g in test.groupby(dimension): rows.append({"dimension":dimension,"group":str(name),"scaffold_count":int(g.scaffold_id.nunique()),**safe_metrics((g.binary_label=="active").astype(int),g.probability)})
    metrics=pd.DataFrame(rows); metrics.to_csv(output/"scaffold_frequency_metrics.csv",index=False); test.assign(error=((test.probability>=.5)!=(test.binary_label=="active"))).to_csv(output/"scaffold_error_analysis.csv",index=False)
    summary={"test_scaffolds":int(test.scaffold_id.nunique()),"unseen_test_scaffolds":int((~test.seen_in_train).groupby(test.scaffold_id).first().sum()),"scaffold_leakage":int((d.groupby("scaffold_id").dataset_split.nunique()>1).sum())}; write_json(output/"scaffold_generalization_summary.json",summary); return metrics
