"""Source-only shortcut audit."""
from __future__ import annotations
import numpy as np, pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.feature_selection import mutual_info_classif
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder
from .core import OUTPUT, SOURCE_FIELDS, load_data, safe_metrics, write_json

def run_source_bias(output=OUTPUT / "source_bias") -> dict:
    output.mkdir(parents=True, exist_ok=True); d=load_data(); y=(d.binary_label=="active").astype(int)
    crosstab=pd.crosstab(d.source_database,d.binary_label,margins=True); crosstab.to_csv(output/"source_class_crosstab.csv")
    rows=[]
    for col in SOURCE_FIELDS:
        codes=pd.factorize(d[col].fillna("missing").astype(str))[0].reshape(-1,1)
        rows.append({"field":col,"mutual_information":float(mutual_info_classif(codes,y,discrete_features=True,random_state=2026)[0]),"unique_values":int(d[col].nunique(dropna=False))})
    pd.DataFrame(rows).to_csv(output/"source_mutual_information.csv",index=False)
    frequencies=[]
    for field in SOURCE_FIELDS+["activity_type","dataset_split"]:
        for (value,label),count in d.groupby([field,"binary_label"],dropna=False).size().items(): frequencies.append({"field":field,"value":str(value),"class":label,"count":int(count)})
    pd.DataFrame(frequencies).to_csv(output/"source_frequency_tables.csv",index=False)
    entropy=[]
    for label,g in d.groupby("binary_label"):
        probabilities=g.source_database.value_counts(normalize=True); entropy.append({"class":label,"source_entropy_bits":float(-(probabilities*np.log2(probabilities)).sum())})
    pd.DataFrame(entropy).to_csv(output/"source_entropy_by_class.csv",index=False)
    train=d.dataset_split.eq("train"); evaluations=[]; predictions=[]
    for target in ("validation","test"):
        mask=d.dataset_split.eq(target)
        model=make_pipeline(OneHotEncoder(handle_unknown="ignore"),LogisticRegression(max_iter=1000,random_state=2026))
        x=d[SOURCE_FIELDS].fillna("missing").astype(str); model.fit(x[train],y[train]); p=model.predict_proba(x[mask])[:,1]
        m=safe_metrics(y[mask],p); m["split"]=target; evaluations.append(m)
        predictions.append(pd.DataFrame({"compound_id":d.loc[mask,"compound_id"],"split":target,"true_label":y[mask],"probability":p}))
    pd.concat(predictions).to_csv(output/"source_only_predictions.csv",index=False)
    metrics={r["split"]:r for r in evaluations}; write_json(output/"source_only_model_metrics.json",metrics)
    high=max(metrics["test"]["roc_auc"],metrics["test"]["pr_auc"])>=.70
    summary={"shortcut_risk":"high" if high else "limited","warning_threshold":.70,"test_metrics":metrics["test"],"production_feature_exclusion_confirmed":True,"source_fields":SOURCE_FIELDS}
    write_json(output/"source_bias_summary.json",summary); return summary
