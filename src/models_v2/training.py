"""Reproducible vector-model training and evaluation."""
from __future__ import annotations
import json,time
from pathlib import Path
import joblib,numpy as np,pandas as pd
from sklearn.metrics import average_precision_score,balanced_accuracy_score,matthews_corrcoef,roc_auc_score
from .backends import create_model
from .representations import fingerprints
from .reproducibility import experiment_metadata
from .schema import validate_dataset
def train_vector_model(dataset_path,output_dir,model_name="random_forest",seed=2026):
 data=pd.read_csv(dataset_path); validate_dataset(data); x=fingerprints(data.canonical_isomeric_smiles); y=data.binary_label.eq("active").astype(int).to_numpy(); train=data.dataset_split.eq("train").to_numpy(); test=data.dataset_split.eq("test").to_numpy(); model=create_model(model_name,seed); start=time.perf_counter(); model.fit(x[train],y[train],sample_weight=data.loc[train,"record_weight"]); duration=time.perf_counter()-start; p=model.predict_proba(x[test])[:,1]; pred=p>=.5; metrics={"roc_auc":float(roc_auc_score(y[test],p)),"pr_auc":float(average_precision_score(y[test],p)),"mcc":float(matthews_corrcoef(y[test],pred)),"balanced_accuracy":float(balanced_accuracy_score(y[test],pred))}; output=Path(output_dir); output.mkdir(parents=True,exist_ok=True); joblib.dump(model,output/f"{model_name}.joblib"); meta=experiment_metadata(dataset_path,seed,model.get_params(),duration); (output/f"{model_name}_metadata.json").write_text(json.dumps({**meta,"metrics":metrics},indent=2)+"\n"); return metrics
