#!/usr/bin/env python3
from __future__ import annotations
import argparse,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.features.feature_matrix import build_feature_artifacts
from src.ml.experiment_runner import run_experiments
from src.explainability.explanation_report import generate_explanations
from src.xai_reports.final_report import generate_xai_report
def main():
 p=argparse.ArgumentParser();
 for flag in ("features-only","train-only","evaluate-only","explain-only","report-only","dashboard-assets-only","overwrite","strict","skip-optional-models","offline"):p.add_argument(f"--{flag}",action="store_true")
 p.add_argument("--seeds",nargs="+",type=int,default=[2026,2027,2028,2029,2030]);p.add_argument("--bootstrap-iterations",type=int,default=500);p.add_argument("--models",nargs="+");p.add_argument("--feature-sets",nargs="+");a=p.parse_args();started=time.monotonic();stages={}
 if not (a.train_only or a.evaluate_only or a.explain_only or a.report_only):
  t=time.monotonic();build_feature_artifacts(ROOT/"data/processed/xai_dataset/hivpr_xai_2000_compounds.csv",ROOT/"data/processed/xai_features",True);stages["features"]=time.monotonic()-t
 if not (a.features_only or a.evaluate_only or a.explain_only or a.report_only):
  t=time.monotonic();run_experiments(ROOT/"data/processed/xai_features",ROOT/"data/processed/xai_models",ROOT/"data/processed/xai_predictions",ROOT/"data/processed/xai_benchmarks",a.seeds,a.models or ("logistic_regression","random_forest","extra_trees","hist_gradient_boosting"),a.feature_sets or ("descriptor","fingerprint","fused"),a.bootstrap_iterations);stages["training_evaluation"]=time.monotonic()-t
 if not (a.features_only or a.train_only or a.evaluate_only or a.report_only):
  t=time.monotonic();generate_explanations(ROOT);stages["explanations"]=time.monotonic()-t
 if not (a.features_only or a.train_only or a.evaluate_only or a.explain_only):
  t=time.monotonic();generate_xai_report(ROOT);stages["reports"]=time.monotonic()-t
 print({"runtime_seconds":time.monotonic()-started,"stages":stages});return 0
if __name__=="__main__":raise SystemExit(main())
