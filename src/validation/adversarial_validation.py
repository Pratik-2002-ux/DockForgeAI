"""Chemical-feature adversarial domain validation."""
from __future__ import annotations
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from .core import OUTPUT, load_data, load_features, safe_metrics

def run_adversarial_validation(output=OUTPUT/"source_bias"):
    d=load_data(); rows=[]; preds=[]; importance=[]
    for family in ("descriptor","fingerprint","fused"):
        x,_,_=load_features(family)
        for name,a,b in (("train_vs_test","train","test"),("train_vs_validation","train","validation")):
            use=d.dataset_split.isin([a,b]).to_numpy(); yy=(d.loc[use,"dataset_split"]==b).astype(int).to_numpy()
            model=make_pipeline(StandardScaler(with_mean=False),LogisticRegression(max_iter=2000,random_state=2026,class_weight="balanced")); cv=StratifiedGroupKFold(n_splits=5,shuffle=True,random_state=2026); p=cross_val_predict(model,x[use],yy,groups=d.loc[use,"scaffold_id"],cv=cv,method="predict_proba")[:,1]; model.fit(x[use],yy)
            m=safe_metrics(yy,p); m.update(comparison=name,feature_set=family); rows.append(m)
        yy=d.source_database.str.lower().eq("chembl").astype(int).to_numpy(); model=make_pipeline(StandardScaler(with_mean=False),LogisticRegression(max_iter=2000,random_state=2026,class_weight="balanced")); cv=StratifiedGroupKFold(n_splits=5,shuffle=True,random_state=2026); p=cross_val_predict(model,x,yy,groups=d.scaffold_id,cv=cv,method="predict_proba")[:,1]; model.fit(x,yy)
        m=safe_metrics(yy,p); m.update(comparison="chembl_vs_dude",feature_set=family); rows.append(m)
        preds.append(pd.DataFrame({"compound_id":d.compound_id,"comparison":"chembl_vs_dude","feature_set":family,"true_domain":yy,"probability":p}))
        coef=np.asarray(model[-1].coef_).ravel(); top=np.argsort(np.abs(coef))[-20:][::-1]
        importance.extend({"feature_set":family,"feature_index":int(i),"importance":float(abs(coef[i]))} for i in top)
    metrics=pd.DataFrame(rows); metrics[metrics.comparison.str.startswith("train")].to_csv(output/"train_test_adversarial_metrics.csv",index=False)
    metrics[metrics.comparison.eq("chembl_vs_dude")].to_csv(output/"source_adversarial_metrics.csv",index=False)
    pd.concat(preds).to_csv(output/"adversarial_predictions.csv",index=False); pd.DataFrame(importance).to_csv(output/"adversarial_feature_importance.csv",index=False)
    return metrics
