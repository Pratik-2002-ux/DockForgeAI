#!/usr/bin/env python3
from __future__ import annotations
import argparse,hashlib,json,platform,subprocess,sys,time
from pathlib import Path
import joblib,numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from src.models_v2.regression.regression_models import models
from src.models_v2.regression.regression_evaluator import metrics
def load_features(root):
 out={"descriptors":pd.read_parquet(root/"descriptors.parquet").to_numpy()}
 for k in ["morgan_r2","morgan_r3","maccs","fused_r2","fused_r3"]:out[k]=np.load(root/(k+".npz"))["X"]
 return out
def main():
 p=argparse.ArgumentParser();p.add_argument("--seeds",nargs="+",type=int,default=[2026,2027,2028,2029,2030]);p.add_argument("--bootstrap-iterations",type=int,default=1000);p.add_argument("--overwrite",action="store_true");a=p.parse_args();d=pd.read_csv(ROOT/"data/model_v2/regression/dataset/ki_regression_dataset.csv");features=load_features(ROOT/"data/model_v2/regression/features");y=d.median_pki.to_numpy();train=d.dataset_split.eq("train").to_numpy();val=d.dataset_split.eq("validation").to_numpy();test=d.dataset_split.eq("test").to_numpy();rows=[];art=ROOT/"data/model_v2/regression/models";ev=ROOT/"data/model_v2/regression/evaluation";art.mkdir(parents=True,exist_ok=True);ev.mkdir(parents=True,exist_ok=True)
 for fname,X in features.items():
  for name,model in models(a.seeds[0]).items():
   start=time.perf_counter();model.fit(X[train],y[train]);vp=model.predict(X[val]);tp=model.predict(X[test]);rows.append({"model":name,"feature_set":fname,"seed":a.seeds[0],"runtime_seconds":time.perf_counter()-start,**{f"validation_{k}":v for k,v in metrics(y[val],vp).items()},**{f"test_{k}":v for k,v in metrics(y[test],tp).items()}})
 comp=pd.DataFrame(rows).sort_values(["validation_rmse","validation_spearman","validation_mae"],ascending=[True,False,True]);comp.to_csv(ev/"model_comparison.csv",index=False);best=comp.iloc[0];X=features[best.feature_set]
 calibration_model=models(a.seeds[0])[best.model];calibration_model.fit(X[train],y[train]);conformal_q=float(np.quantile(np.abs(y[val]-calibration_model.predict(X[val])),.95))
 ensemble=[];best_model=None
 for seed in a.seeds:
  best_model=models(seed)[best.model];best_model.fit(X[train|val],y[train|val]);ensemble.append(best_model.predict(X[test]))
 joblib.dump(best_model,art/"best_ki_regressor.joblib");ensemble=np.asarray(ensemble);pred=ensemble.mean(axis=0);uncertainty=ensemble.std(axis=0)
 fp=features["morgan_r3"].astype(bool);train_fp=fp[train|val];test_fp=fp[test];intersection=test_fp.astype(int)@train_fp.astype(int).T;union=test_fp.sum(1)[:,None]+train_fp.sum(1)[None,:]-intersection;similarity=np.divide(intersection,union,out=np.zeros_like(intersection,dtype=float),where=union>0);max_similarity=similarity.max(1);top5=np.sort(similarity,axis=1)[:,-5:].mean(1);domain=max_similarity>=.30
 pd.DataFrame({"compound_id":d.loc[test,"compound_id"],"observed_pki":y[test],"predicted_pki":pred,"ensemble_std":uncertainty,"interval_lower":pred-conformal_q,"interval_upper":pred+conformal_q,"maximum_training_tanimoto":max_similarity,"mean_top5_tanimoto":top5,"applicability_domain":np.where(domain,"in_domain","out_of_domain")}).to_csv(ev/"locked_scaffold_test_predictions.csv",index=False)
 rng=np.random.default_rng(2026);boot=[]
 for _ in range(a.bootstrap_iterations):
  idx=rng.integers(0,len(pred),len(pred));boot.append(metrics(y[test][idx],pred[idx])["rmse"])
 result={"selection_uses_external":False,"best_model":best.model,"best_feature_set":best.feature_set,"validation_metrics":{k.removeprefix("validation_"):v for k,v in best.items() if k.startswith("validation_")},"locked_scaffold_test_metrics":metrics(y[test],pred),"test_rmse_ci95":[float(np.quantile(boot,.025)),float(np.quantile(boot,.975))],"uncertainty":{"ensemble_seeds":a.seeds,"mean_prediction_std":float(uncertainty.mean()),"conformal_half_width_95":conformal_q},"applicability_domain":{"in_domain_count":int(domain.sum()),"out_of_domain_count":int((~domain).sum()),"threshold_max_tanimoto":.30},"model_path":"data/model_v2/regression/models/best_ki_regressor.joblib"};(ev/"best_model.json").write_text(json.dumps(result,indent=2,default=float)+"\n")
 model_path=art/"best_ki_regressor.joblib";dataset_path=ROOT/"data/model_v2/regression/dataset/ki_regression_dataset.csv";git_hash=subprocess.run(["git","rev-parse","HEAD"],cwd=ROOT,text=True,capture_output=True,check=True).stdout.strip();training_manifest={"model_id":f"ki_rf_morgan_r3_{git_hash[:8]}","endpoint":"exact Ki converted to pKi","model_family":best.model,"feature_set":best.feature_set,"seeds":a.seeds,"dataset_sha256":hashlib.sha256(dataset_path.read_bytes()).hexdigest(),"model_sha256":hashlib.sha256(model_path.read_bytes()).hexdigest(),"git_commit":git_hash,"python":platform.python_version(),"platform":platform.platform(),"external_benchmark_used_for_selection":False,"validation_selection_metric":"scaffold validation RMSE","training_rows":int(train.sum()),"validation_rows":int(val.sum()),"test_rows":int(test.sum())};(art/"training_manifest.json").write_text(json.dumps(training_manifest,indent=2)+"\n");model_card={"model_id":training_manifest["model_id"],"intended_use":"Research-only pKi prediction for structures within the represented chemical domain.","not_intended_for":["clinical decisions","therapeutic claims","molecule optimization","binary activity classification"],"dataset_limitations":"Small, BindingDB-only, 192-compound exact-Ki dataset; scaffold test n=20.","metrics":result,"explanations":"Describe model behaviour and are not biological causality."};(art/"model_card.json").write_text(json.dumps(model_card,indent=2,default=float)+"\n");print(json.dumps(result,indent=2,default=float));return 0
if __name__=="__main__":raise SystemExit(main())
