#!/usr/bin/env python3
"""Generate conservative feature importance for the frozen Ki model."""
import json
from pathlib import Path
import joblib,numpy as np,pandas as pd,shap
ROOT=Path(__file__).resolve().parents[2]
if __name__=="__main__":
 m=joblib.load(ROOT/"data/model_v2/regression/models/best_ki_regressor.joblib");out=ROOT/"data/model_v2/regression/explanations";out.mkdir(parents=True,exist_ok=True);d=pd.read_csv(ROOT/"data/model_v2/regression/dataset/ki_regression_dataset.csv");X=np.load(ROOT/"data/model_v2/regression/features/morgan_r3.npz")["X"];test=d.dataset_split.eq("test").to_numpy();train=~test;values=getattr(m,"feature_importances_",np.array([]));pd.DataFrame({"feature_index":range(len(values)),"importance":values}).sort_values("importance",ascending=False).to_csv(out/"feature_importance.csv",index=False)
 sv=shap.TreeExplainer(m).shap_values(X[test]);np.savez_compressed(out/"shap_values.npz",values=np.asarray(sv),compound_ids=d.loc[test,"compound_id"].to_numpy())
 rng=np.random.default_rng(2026);baseline=np.sqrt(np.mean((d.loc[test,"median_pki"].to_numpy()-m.predict(X[test]))**2));perm_rows=[]
 for feature in np.argsort(values)[-50:]:
  scores=[]
  for _ in range(10):
   changed=X[test].copy();changed[:,feature]=rng.permutation(changed[:,feature]);scores.append(np.sqrt(np.mean((d.loc[test,"median_pki"].to_numpy()-m.predict(changed))**2))-baseline)
  perm_rows.append({"feature_index":int(feature),"importance_mean":float(np.mean(scores)),"importance_std":float(np.std(scores))})
 pd.DataFrame(perm_rows).sort_values("importance_mean",ascending=False).to_csv(out/"permutation_importance.csv",index=False)
 xb=X.astype(bool);inter=xb[test].astype(int)@xb[train].astype(int).T;union=xb[test].sum(1)[:,None]+xb[train].sum(1)[None,:]-inter;sim=np.divide(inter,union,out=np.zeros_like(inter,dtype=float),where=union>0);train_ids=d.loc[train,"compound_id"].to_numpy();rows=[]
 for compound,sims in zip(d.loc[test,"compound_id"],sim):
  for rank,index in enumerate(np.argsort(sims)[-5:][::-1],1):rows.append({"compound_id":compound,"neighbor_rank":rank,"neighbor_compound_id":train_ids[index],"tanimoto":sims[index]})
 pd.DataFrame(rows).to_csv(out/"nearest_neighbors.csv",index=False);(out/"explanation_manifest.json").write_text(json.dumps({"describes_model_behavior_only":True,"biological_causality":False,"hash_collisions_possible":True,"artifacts":["feature_importance.csv","permutation_importance.csv","shap_values.npz","nearest_neighbors.csv"]},indent=2)+"\n")
