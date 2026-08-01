"""Permutation and Y-scrambling controls."""
from __future__ import annotations
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from .core import OUTPUT, best_predictions, load_data, load_features, safe_metrics, write_json

def empirical_p_value(scores, observed): return float((np.sum(np.asarray(scores)>=observed)+1)/(len(scores)+1))

def run_permutations(n=100, output=OUTPUT/"permutation"):
    x,y,index=load_features("descriptor"); d=load_data(); tr=index.dataset_split.eq("train").to_numpy(); va=index.dataset_split.eq("validation").to_numpy(); te=index.dataset_split.eq("test").to_numpy(); rng=np.random.default_rng(2026)
    rows=[]
    for i in range(n):
        yp=rng.permutation(y[tr]); model=LogisticRegression(max_iter=1000,random_state=2026).fit(x[tr],yp)
        for split,mask in (("validation",va),("test",te)):
            m=safe_metrics(y[mask],model.predict_proba(x[mask])[:,1]); rows.append({"repetition":i,"split":split,**m})
    frame=pd.DataFrame(rows); frame.to_csv(output/"label_permutation_metrics.csv",index=False)
    for kind,group in (("scaffold",d.scaffold_id),("source",d.source_database)):
        block=y.copy()
        for _,ids in d.groupby(group).groups.items(): block[list(ids)]=rng.permutation(block[list(ids)])
        model=LogisticRegression(max_iter=1000,random_state=2026).fit(x[tr],block[tr]); m=safe_metrics(y[te],model.predict_proba(x[te])[:,1])
        pd.DataFrame([{**m,"repetition":0}]).to_csv(output/f"{kind}_block_permutation_metrics.csv",index=False)
    observed=safe_metrics((best_predictions().true_label=="active").astype(int),best_predictions().predicted_probability)["pr_auc"]
    scores=frame.loc[frame.split.eq("test"),"pr_auc"].to_numpy(); summary={"permutations":n,"observed_test_pr_auc":observed,"permuted_test_pr_auc_mean":float(scores.mean()),"permuted_test_pr_auc_std":float(scores.std(ddof=1)),"empirical_p_value":empirical_p_value(scores,observed),"full_y_scrambling":"diagnostic training-label permutation used; locked labels retained for scoring"}
    write_json(output/"permutation_summary.json",summary); return summary
