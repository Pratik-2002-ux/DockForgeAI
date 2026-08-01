from __future__ import annotations
from sklearn.ensemble import ExtraTreesClassifier,HistGradientBoostingClassifier,RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
def create_model(name:str,seed:int,params:dict|None=None):
 p=params or {}
 if name=="logistic_regression": return make_pipeline(StandardScaler(with_mean=False),LogisticRegression(C=p.get("C",1.0),max_iter=2000,class_weight="balanced",random_state=seed))
 if name=="random_forest": return RandomForestClassifier(n_estimators=p.get("n_estimators",200),max_depth=p.get("max_depth",None),min_samples_leaf=p.get("min_samples_leaf",2),max_features="sqrt",class_weight="balanced",n_jobs=-1,random_state=seed)
 if name=="extra_trees": return ExtraTreesClassifier(n_estimators=p.get("n_estimators",200),max_depth=p.get("max_depth",None),min_samples_leaf=p.get("min_samples_leaf",2),class_weight="balanced",n_jobs=-1,random_state=seed)
 if name=="hist_gradient_boosting": return HistGradientBoostingClassifier(learning_rate=p.get("learning_rate",.1),max_iter=p.get("max_iter",100),max_leaf_nodes=31,l2_regularization=1,random_state=seed)
 raise ValueError(f"Unknown model: {name}")
