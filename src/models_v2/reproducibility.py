"""Experiment metadata capture."""
from __future__ import annotations
import hashlib,json,platform,subprocess,sys,time
from pathlib import Path
import importlib.metadata
def checksum(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def git_hash():
 try: return subprocess.check_output(["git","rev-parse","HEAD"],text=True).strip()
 except Exception: return "unavailable"
def package_versions(names=("numpy","pandas","scikit-learn","rdkit","shap")):
 result={}
 for name in names:
  try: result[name]=importlib.metadata.version(name)
  except importlib.metadata.PackageNotFoundError: result[name]="not_installed"
 return result
def experiment_metadata(dataset_path,seed,hyperparameters,duration,device="cpu"):
 return {"random_seed":seed,"package_versions":package_versions(),"dataset_checksum":checksum(dataset_path),"git_hash":git_hash(),"hyperparameters":hyperparameters,"training_duration_seconds":float(duration),"device":device,"cpu":platform.processor() or platform.machine(),"platform":platform.platform(),"python":sys.version.split()[0]}
