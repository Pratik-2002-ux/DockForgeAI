"""Non-parametric chemical-similarity baselines."""
from __future__ import annotations
import numpy as np, pandas as pd
from sklearn.metrics import pairwise_distances
from .core import OUTPUT, jaccard_similarity, load_features, safe_metrics, write_json

def run_nearest_neighbor(output=OUTPUT/"baselines"):
    fp,y,index=load_features("fingerprint"); desc,_,_=load_features("descriptor"); tr=index.dataset_split.eq("train").to_numpy(); rows=[]; metrics=[]
    for split in ("validation","test"):
        mask=index.dataset_split.eq(split).to_numpy(); sim=jaccard_similarity(fp[mask],fp[tr]); order=np.argsort(sim,axis=1)[:,::-1]; train_y=y[tr]
        for k in (1,3,5):
            ids=order[:,:k]; probs=train_y[ids].mean(axis=1); consistency=np.maximum(probs,1-probs)
            metrics.append({"baseline":f"tanimoto_{k}nn","split":split,**safe_metrics(y[mask],probs)})
            rows.append(pd.DataFrame({"compound_id":index.loc[mask,"compound_id"].to_numpy(),"split":split,"baseline":f"tanimoto_{k}nn","probability":probs,"maximum_similarity":sim.max(1),"mean_top_k_similarity":np.take_along_axis(sim,ids,axis=1).mean(1),"neighbour_label_consistency":consistency}))
        ids=order[:,:5]; weights=np.take_along_axis(sim,ids,axis=1); probs=(train_y[ids]*weights).sum(1)/np.maximum(weights.sum(1),1e-12)
        metrics.append({"baseline":"similarity_weighted_5nn","split":split,**safe_metrics(y[mask],probs)})
        dist=pairwise_distances(desc[mask],desc[tr]); ids=np.argmin(dist,axis=1); probs=train_y[ids].astype(float); metrics.append({"baseline":"descriptor_1nn","split":split,**safe_metrics(y[mask],probs)})
    pred=pd.concat(rows,ignore_index=True); met=pd.DataFrame(metrics); pred.to_csv(output/"nearest_neighbor_predictions.csv",index=False); met.to_csv(output/"nearest_neighbor_metrics.csv",index=False)
    best=met[met.split.eq("test")].sort_values("pr_auc",ascending=False).iloc[0].to_dict(); summary={"best_test_baseline":best,"local_similarity_warning":bool(best["pr_auc"]>=.95),"interpretation":"Near-best nearest-neighbour performance indicates local similarity may dominate."}; write_json(output/"nearest_neighbor_summary.json",summary); return met
