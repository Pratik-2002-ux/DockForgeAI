"""Explicit model backend availability and construction."""
from __future__ import annotations
import importlib.util
from sklearn.ensemble import RandomForestClassifier
BACKENDS={"random_forest":("sklearn",True),"xgboost":("xgboost",False),"lightgbm":("lightgbm",False),"catboost":("catboost",False),"chemprop_mpnn":("chemprop",False),"gin":("torch_geometric",False),"gine":("torch_geometric",False),"attentivefp":("torch_geometric",False),"graphgps":("torch_geometric",False),"gin_transformer":("torch_geometric",False)}
class BackendUnavailable(RuntimeError): pass
def availability(): return {name:(required or importlib.util.find_spec(module) is not None) for name,(module,required) in BACKENDS.items()}
def create_model(name,seed=2026):
 if name=="random_forest": return RandomForestClassifier(n_estimators=300,min_samples_leaf=2,class_weight="balanced",random_state=seed,n_jobs=-1)
 if not availability().get(name,False): raise BackendUnavailable(f"{name} requires optional backend {BACKENDS[name][0]}; no substitute was used")
 if name=="xgboost":
  from xgboost import XGBClassifier; return XGBClassifier(n_estimators=300,random_state=seed,n_jobs=-1)
 if name=="lightgbm":
  from lightgbm import LGBMClassifier; return LGBMClassifier(n_estimators=300,random_state=seed)
 if name=="catboost":
  from catboost import CatBoostClassifier; return CatBoostClassifier(iterations=300,random_seed=seed,verbose=False)
 raise BackendUnavailable(f"{name} adapter requires its graph training entrypoint")
