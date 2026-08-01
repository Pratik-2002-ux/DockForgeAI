"""Source-balanced validation models without modifying locked artifacts."""
from __future__ import annotations
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from .core import OUTPUT, load_data, load_features, safe_metrics, write_json

def source_weights(data, method="joint", clip=(.1,10.)):
    key=data.source_database if method=="source" else data.source_database.astype(str)+"|"+data.binary_label.astype(str); w=len(data)/(key.nunique()*key.map(key.value_counts())); return w.clip(*clip).to_numpy(float)
def run_reweighting(output=OUTPUT/"reweighting"):
    d=load_data(); x,y,index=load_features("descriptor"); tr=index.dataset_split.eq("train").to_numpy(); te=index.dataset_split.eq("test").to_numpy(); rows=[]; wd=[]
    for method in ("unweighted","inverse_source_frequency","class_source_joint"):
        w=None if method=="unweighted" else source_weights(d[tr],"source" if method.startswith("inverse") else "joint"); model=LogisticRegression(max_iter=1000,random_state=2026).fit(x[tr],y[tr],sample_weight=w); rows.append({"method":method,**safe_metrics(y[te],model.predict_proba(x[te])[:,1])})
        if w is not None: wd.extend({"method":method,"weight":float(v)} for v in w)
    frame=pd.DataFrame(rows); frame.to_csv(output/"source_weighted_metrics.csv",index=False); pd.DataFrame(wd).to_csv(output/"weight_distribution.csv",index=False); summary={"models_are_validation_only":True,"locked_artifacts_overwritten":False,"metrics":rows}; write_json(output/"reweighting_summary.json",summary); return frame
