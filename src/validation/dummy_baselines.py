"""Locked-test dummy baselines."""
from __future__ import annotations
import numpy as np, pandas as pd
from sklearn.dummy import DummyClassifier
from .core import OUTPUT, load_features, safe_metrics

def run_dummy_baselines(output=OUTPUT/"baselines"):
    x,y,index=load_features("descriptor"); tr=index.dataset_split.eq("train").to_numpy(); te=index.dataset_split.eq("test").to_numpy(); rows=[]
    for strategy in ("most_frequent","prior","stratified","uniform"):
        model=DummyClassifier(strategy=strategy,random_state=2026).fit(x[tr],y[tr]); p=model.predict_proba(x[te]); p=p[:,list(model.classes_).index(1)] if 1 in model.classes_ else np.zeros(te.sum()); rows.append({"baseline":strategy,**safe_metrics(y[te],p)})
    rows.append({"baseline":"constant_0.5",**safe_metrics(y[te],np.full(te.sum(),.5))}); frame=pd.DataFrame(rows); frame.to_csv(output/"dummy_baselines.csv",index=False); return frame
