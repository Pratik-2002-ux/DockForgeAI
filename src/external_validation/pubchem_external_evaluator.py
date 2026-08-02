"""Locked-model evaluation of the curated PubChem dataset."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.spatial.distance import jensenshannon
from scipy.stats import ks_2samp, wasserstein_distance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold,cross_val_predict
from src.ml.metrics import classification_metrics
from src.features.descriptor_features import descriptor_frame
from src.validation.external_dataset_validator import validate_external_dataset

def bootstrap_metrics(y,p,iterations=1000,seed=2026):
    rng=np.random.default_rng(seed); rows=[]
    for _ in range(iterations):
        ids=rng.integers(0,len(y),len(y)); yy=np.asarray(y)[ids]
        if len(np.unique(yy))<2: continue
        m=classification_metrics(yy,np.asarray(p)[ids]); rows.append({k:m[k] for k in ("roc_auc","pr_auc","mcc")})
    return {k:{"lower_95":float(pd.DataFrame(rows)[k].quantile(.025)),"upper_95":float(pd.DataFrame(rows)[k].quantile(.975)),"valid_iterations":len(rows)} for k in ("roc_auc","pr_auc","mcc")} if rows else {}
def group_metrics(frame,iterations):
    y=frame.true_label.eq("active").astype(int).to_numpy(); p=frame.predicted_probability.to_numpy(); result={"sample_count":len(frame),"active_count":int(y.sum()),"inactive_count":int((1-y).sum()),"unstable":len(frame)<30 or len(np.unique(y))<2 or min(y.sum(),(1-y).sum())<5}
    if len(np.unique(y))==2: result.update(classification_metrics(y,p)); result["bootstrap_95_ci"]=bootstrap_metrics(y,p,iterations) if not result["unstable"] else {}
    return result
def evaluate_locked_model(dataset_path:Path,output_dir:Path,bootstrap_iterations=1000,overwrite=False):
    output=Path(output_dir); output.mkdir(parents=True,exist_ok=True); data=pd.read_csv(dataset_path); inference=output/"locked_inference"; validate_external_dataset(dataset_path,"compound_id","canonical_isomeric_smiles",inference,"true_label",strict=True,allow_scaffold_overlap=True,overwrite=overwrite); pred=pd.read_csv(inference/"external_predictions.csv"); merged=data.merge(pred[["compound_id","predicted_probability","predicted_class"]],on="compound_id",validate="one_to_one"); merged.to_csv(output/"pubchem_external_predictions.csv",index=False)
    evaluations={"all_exact_independent":group_metrics(merged,bootstrap_iterations)}
    groups={"scaffold_overlapping":merged.scaffold_training_overlap,"novel_scaffold":~merged.scaffold_training_overlap,"in_domain":merged.applicability_domain_status.eq("in_domain"),"out_of_domain":merged.applicability_domain_status.eq("out_of_domain")}
    for name,mask in groups.items(): evaluations[name]=group_metrics(merged[mask],bootstrap_iterations)
    per={str(int(aid)):group_metrics(g,bootstrap_iterations) for aid,g in merged.groupby("pubchem_aid")}; (output/"external_validation_metrics.json").write_text(json.dumps({"locked_model":True,"threshold":0.5,"threshold_optimized":False,"model_retrained":False,"groups":evaluations,"per_assay":per},indent=2,default=lambda x:x.item() if hasattr(x,"item") else str(x))+"\n")
    pd.DataFrame([{"group":k,**{x:y for x,y in v.items() if not isinstance(y,dict)}} for k,v in evaluations.items()]).to_csv(output/"group_metrics.csv",index=False); pd.DataFrame([{"pubchem_aid":k,**{x:y for x,y in v.items() if not isinstance(y,dict)}} for k,v in per.items()]).to_csv(output/"per_assay_metrics.csv",index=False)
    training=pd.read_csv("data/processed/xai_dataset/hivpr_xai_2000_compounds.csv"); columns=["molecular_weight","logp","tpsa","hbd","hba","rotatable_bonds","formal_charge","heavy_atom_count"]; ext_desc=descriptor_frame(merged); shifts=[]
    for c in columns:
        x=training[c].dropna().to_numpy(); y=ext_desc[c].dropna().to_numpy(); shifts.append({"descriptor":c,"ks_statistic":ks_2samp(x,y).statistic,"wasserstein_distance":wasserstein_distance(x,y)})
    pd.DataFrame(shifts).to_csv(output/"distribution_shift.csv",index=False); combined=np.vstack([training[columns].to_numpy(),ext_desc[columns].to_numpy()]); domain=np.r_[np.zeros(len(training),int),np.ones(len(ext_desc),int)]; cv=StratifiedKFold(5,shuffle=True,random_state=2026); probability=cross_val_predict(LogisticRegression(max_iter=2000),combined,domain,cv=cv,method="predict_proba")[:,1]; adversarial=float(roc_auc_score(domain,probability)); summary={"adversarial_train_external_roc_auc":adversarial,"distribution_shift_warning":adversarial>=.70,"applicability_domain_coverage":float(merged.applicability_domain_status.eq("in_domain").mean()),"exact_overlap_removed":True,"external_validation_status":"completed_with_caution","locked_model_preserved":True}; (output/"external_validation_summary.json").write_text(json.dumps(summary,indent=2)+"\n"); return summary
