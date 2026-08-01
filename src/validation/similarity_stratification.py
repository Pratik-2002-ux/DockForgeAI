"""Performance stratified by training-set similarity."""
from __future__ import annotations
import numpy as np, pandas as pd
from .core import OUTPUT, best_predictions, jaccard_similarity, load_features, safe_metrics

def run_similarity_stratification(output=OUTPUT/"similarity"):
    x,y,index=load_features("fingerprint"); tr=index.dataset_split.eq("train").to_numpy(); te=index.dataset_split.eq("test").to_numpy(); sim=jaccard_similarity(x[te],x[tr]).max(1); pred=best_predictions(); bins=[0,.2,.4,.6,.8,1.0000001]; labels=["[0.0,0.2)","[0.2,0.4)","[0.4,0.6)","[0.6,0.8)","[0.8,1.0]"]; groups=pd.cut(sim,bins,labels=labels,include_lowest=True,right=False); dist=pd.DataFrame({"compound_id":index.loc[te,"compound_id"].to_numpy(),"maximum_similarity":sim,"similarity_bin":groups.astype(str),"true_label":y[te],"probability":pred.predicted_probability})
    dist.to_csv(output/"similarity_distribution.csv",index=False); rows=[]
    for name,g in dist.groupby("similarity_bin",observed=True): rows.append({"similarity_bin":name,**safe_metrics(g.true_label,g.probability)})
    pd.DataFrame(rows).to_csv(output/"similarity_stratified_metrics.csv",index=False); dist.assign(error=((dist.probability>=.5)!=dist.true_label)).to_csv(output/"similarity_error_analysis.csv",index=False); return pd.DataFrame(rows)
