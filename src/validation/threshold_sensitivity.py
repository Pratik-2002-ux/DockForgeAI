"""Cached-record activity-threshold sensitivity."""
from __future__ import annotations
import pandas as pd
from sklearn.linear_model import LogisticRegression
from .core import OUTPUT, load_data, load_features, safe_metrics, write_json

def run_threshold_sensitivity(output=OUTPUT/"thresholds"):
    d=load_data(); x,_,index=load_features("descriptor"); rows=[]; counts=[]
    for threshold in (5.5,6.,6.5,7.):
        y=((d.source_database.eq("ChEMBL") & d.pchembl_value.ge(threshold)) | ~d.source_database.eq("ChEMBL")).astype(int).to_numpy(copy=True); y[~d.source_database.eq("ChEMBL").to_numpy()]=0; eligible=(~d.source_database.eq("ChEMBL") | d.pchembl_value.notna()).to_numpy(); tr=index.dataset_split.eq("train").to_numpy()&eligible; va=index.dataset_split.eq("validation").to_numpy()&eligible; te=index.dataset_split.eq("test").to_numpy()&eligible
        model=LogisticRegression(max_iter=1000,random_state=2026).fit(x[tr],y[tr])
        for split,mask in (("validation",va),("test",te)): rows.append({"threshold":threshold,"split":split,**safe_metrics(y[mask],model.predict_proba(x[mask])[:,1])})
        counts.append({"threshold":threshold,"active":int(y[eligible].sum()),"decoy":int((y[eligible]==0).sum()),"ambiguous_excluded":int((d.source_database.eq("ChEMBL")&d.pchembl_value.lt(threshold)).sum())})
    pd.DataFrame(rows).to_csv(output/"threshold_sensitivity_metrics.csv",index=False); pd.DataFrame(counts).to_csv(output/"threshold_class_counts.csv",index=False); summary={"thresholds_not_selected_by_score":True,"cached_records_only":True,"metrics":rows}; write_json(output/"threshold_summary.json",summary); return pd.DataFrame(rows)
