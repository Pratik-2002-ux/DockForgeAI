"""Run repeated-seed scaffold-safe model experiments."""
from __future__ import annotations
import hashlib,json,platform,subprocess,time
from pathlib import Path
import joblib,numpy as np,pandas as pd
from scipy import sparse
from src.ml.calibration import choose_calibration
from src.ml.dataset_loader import load_feature_set,split_indices
from src.ml.metrics import classification_metrics
from src.ml.model_registry import write_registry
from src.ml.trainer import train_model,_usable

MODELS=("logistic_regression","random_forest","extra_trees","hist_gradient_boosting")
FEATURES=("descriptor","fingerprint","fused")
def run_experiments(feature_dir:Path,model_dir:Path,prediction_dir:Path,benchmark_dir:Path,seeds:list[int],models=MODELS,feature_sets=FEATURES,bootstrap_iterations=500):
 for path in (model_dir,prediction_dir,benchmark_dir): path.mkdir(parents=True,exist_ok=True)
 registry=[]; metric_rows=[]; predictions=[]
 for feature in feature_sets:
  x,y,index=load_feature_set(feature_dir,feature); splits=split_indices(index)
  for seed in seeds:
   for name in models:
    model_id=f"{name}__{feature}__seed{seed}"; artifact=model_dir/f"{model_id}.joblib"
    model,val_metrics,runtime=train_model(name,seed,x,y,splits,artifact)
    val=splits["validation"]; test=splits["test"]; val_raw=model.predict_proba(_usable(x[val],name))[:,1]; method,cal,val_cal_metrics=choose_calibration(y[val],val_raw); test_raw=model.predict_proba(_usable(x[test],name))[:,1]; test_cal=cal.predict(test_raw); joblib.dump(cal,model_dir/f"{model_id}__calibrator.joblib")
    for split,indices,raw,calibrated in (("validation",val,val_raw,cal.predict(val_raw)),("test",test,test_raw,test_cal)):
     metrics=classification_metrics(y[indices],calibrated); metric_rows.append({"model_id":model_id,"model_name":name,"feature_set":feature,"seed":seed,"dataset_split":split,"calibration_method":method,**metrics})
     for pos,idx in enumerate(indices): predictions.append({"compound_id":index.iloc[idx].compound_id,"standardized_compound_id":index.iloc[idx].standardized_compound_id,"dataset_split":split,"true_label":index.iloc[idx].binary_label,"raw_probability":float(raw[pos]),"predicted_probability":float(calibrated[pos]),"predicted_class":"active" if calibrated[pos]>=.5 else "inactive_or_decoy","model_id":model_id,"feature_set":feature,"seed":seed,"scaffold_id":index.iloc[idx].scaffold_id})
    parameters={k:v for k,v in model.get_params().items() if isinstance(v,(str,int,float,bool,type(None)))}
    registry.append({"model_id":model_id,"model_name":name,"feature_set":feature,"seed":seed,"hyperparameters":parameters,"calibration_method":method,"training_timestamp":pd.Timestamp.utcnow().isoformat(),"training_runtime_seconds":runtime,"artifact_path":artifact.relative_to(model_dir.parent.parent.parent).as_posix(),"metrics_path":"data/processed/xai_benchmarks/model_comparison.csv","dataset_checksum":json.loads((feature_dir/"feature_manifest.json").read_text())["dataset_sha256"],"feature_checksum":hashlib.sha256((feature_dir/"feature_manifest.json").read_bytes()).hexdigest(),"git_commit_hash":subprocess.run(["git","rev-parse","--short","HEAD"],capture_output=True,text=True).stdout.strip()})
 metrics_frame=pd.DataFrame(metric_rows); metrics_frame.to_csv(benchmark_dir/"model_comparison.csv",index=False); pred=pd.DataFrame(predictions); pred.to_csv(prediction_dir/"all_model_predictions.csv",index=False); pred[pred.dataset_split=="validation"].to_csv(prediction_dir/"validation_predictions.csv",index=False); pred[pred.dataset_split=="test"].to_csv(prediction_dir/"test_predictions.csv",index=False)
 validation=metrics_frame[metrics_frame.dataset_split=="validation"]; best_row=validation.sort_values(["pr_auc","roc_auc","mcc"],ascending=False).iloc[0]; best_id=best_row.model_id; best_test=metrics_frame[(metrics_frame.model_id==best_id)&(metrics_frame.dataset_split=="test")].iloc[0]
 best={"model_id":best_id,"selection_metric":"validation_pr_auc","validation_metrics":best_row.to_dict(),"locked_test_metrics":best_test.to_dict()}; (benchmark_dir/"best_model.json").write_text(json.dumps(best,indent=2,default=str)); write_registry(model_dir/"model_registry.json",registry)
 summary=validation.groupby(["model_name","feature_set"])[["pr_auc","roc_auc","mcc","brier_score"]].agg(["mean","std"]).reset_index(); summary.to_csv(benchmark_dir/"seed_stability.csv",index=False)
 family=metrics_frame[metrics_frame.dataset_split=="test"].groupby("feature_set")[["pr_auc","roc_auc","mcc","brier_score"]].agg(["mean","std"]).reset_index();family.to_csv(benchmark_dir/"feature_family_comparison.csv",index=False);family.to_csv(benchmark_dir/"ablation_results.csv",index=False)
 pd.DataFrame([{"comparison":"scaffold split (primary)","status":"reported"},{"comparison":"random split (secondary)","status":"not used for primary selection; source shortcut caveat"}]).to_csv(benchmark_dir/"random_vs_scaffold_split.csv",index=False)
 test_metrics=metrics_frame[metrics_frame.dataset_split=="test"];test_metrics.groupby("calibration_method")[["brier_score","expected_calibration_error","pr_auc"]].mean().reset_index().to_csv(benchmark_dir/"calibration_summary.csv",index=False)
 confusion={row.model_id:{k:int(row[k]) for k in ("tn","fp","fn","tp")} for _,row in test_metrics.iterrows()};(benchmark_dir/"confusion_matrices.json").write_text(json.dumps(confusion,indent=2))
 comparisons=[];base=test_metrics.groupby("model_name").pr_auc.mean().sort_values(ascending=False)
 for name,value in base.items(): comparisons.append({"model_name":name,"mean_test_pr_auc":value,"difference_from_best":value-base.iloc[0],"interpretation":"descriptive repeated-seed comparison; source shortcut caveat applies"})
 pd.DataFrame(comparisons).to_csv(benchmark_dir/"model_statistical_comparison.csv",index=False)
 best_predictions=pred[(pred.model_id==best_id)&(pred.dataset_split=="test")];rng=np.random.default_rng(2026);boots=[]
 truth=(best_predictions.true_label=="active").astype(int).to_numpy();prob=best_predictions.predicted_probability.to_numpy()
 for _ in range(bootstrap_iterations):
  indices=rng.integers(0,len(truth),len(truth))
  if len(np.unique(truth[indices]))<2: continue
  metric=classification_metrics(truth[indices],prob[indices]);boots.append({"roc_auc":metric["roc_auc"],"pr_auc":metric["pr_auc"],"mcc":metric["mcc"]})
 pd.DataFrame([{"metric":column,"lower_95":np.quantile(pd.DataFrame(boots)[column],.025),"upper_95":np.quantile(pd.DataFrame(boots)[column],.975),"iterations":len(boots)} for column in ("roc_auc","pr_auc","mcc")]).to_csv(benchmark_dir/"confidence_intervals.csv",index=False)
 (benchmark_dir/"benchmark_manifest.json").write_text(json.dumps({"models":list(models),"feature_sets":list(feature_sets),"seeds":seeds,"primary_split":"scaffold","primary_selection_metric":"validation_pr_auc","test_set_used_for_tuning":False,"bootstrap_iterations":bootstrap_iterations},indent=2)); return best
