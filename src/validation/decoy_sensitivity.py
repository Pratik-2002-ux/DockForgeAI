"""Deterministic sensitivity analysis over the cached DUD-E decoy pool."""
from __future__ import annotations
import numpy as np, pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupShuffleSplit
from sklearn.preprocessing import StandardScaler
from .core import OUTPUT, ROOT, load_data, safe_metrics, write_json

COLUMNS=["molecular_weight","logp","tpsa","hbd","hba","rotatable_bonds","formal_charge","heavy_atom_count"]
def _properties(smiles):
    mol=Chem.MolFromSmiles(smiles)
    if mol is None: return None
    return [Descriptors.MolWt(mol),Crippen.MolLogP(mol),rdMolDescriptors.CalcTPSA(mol),Lipinski.NumHDonors(mol),Lipinski.NumHAcceptors(mol),Lipinski.NumRotatableBonds(mol),Chem.GetFormalCharge(mol),mol.GetNumHeavyAtoms()],MurckoScaffold.MurckoScaffoldSmiles(mol=mol,includeChirality=True)

def run_decoy_sensitivity(seeds=range(2026,2036), output=OUTPUT/"decoy_sensitivity"):
    d=load_data(); active=d[d.binary_label.eq("active")].copy(); active_x=active[COLUMNS].to_numpy(float); active_groups=active.scaffold_id.to_numpy(); pool=[]
    for line in (ROOT/"data/raw/xai/dude_hivpr/decoys_final.ism").read_text().splitlines():
        parts=line.split()
        if len(parts)>=2: pool.append((parts[0],parts[1]))
    rows=[]; props=[]
    for seed in seeds:
        rng=np.random.default_rng(seed); candidate_ids=rng.choice(len(pool),size=min(3000,len(pool)),replace=False); records=[]
        for idx in candidate_ids:
            values=_properties(pool[idx][0])
            if values: records.append((pool[idx][1],*values[0],values[1]))
        candidates=pd.DataFrame(records,columns=["compound_id",*COLUMNS,"scaffold"]); scaler=StandardScaler().fit(pd.concat([active[COLUMNS],candidates[COLUMNS]])); center=scaler.transform(active[COLUMNS]).mean(0); candidates["match_distance"]=np.linalg.norm(scaler.transform(candidates[COLUMNS])-center,axis=1); selected=candidates.nsmallest(min(1000,len(candidates)),"match_distance")
        xx=np.vstack([active_x,selected[COLUMNS].to_numpy(float)]); yy=np.r_[np.ones(len(active),int),np.zeros(len(selected),int)]; groups=np.r_[active_groups,selected.scaffold.to_numpy()]; first=GroupShuffleSplit(n_splits=1,test_size=.30,random_state=seed); train,temp=next(first.split(xx,yy,groups)); second=GroupShuffleSplit(n_splits=1,test_size=.50,random_state=seed); val_rel,test_rel=next(second.split(xx[temp],yy[temp],groups[temp])); test=temp[test_rel]
        model=LogisticRegression(max_iter=1000,random_state=seed,class_weight="balanced").fit(xx[train],yy[train]); m=safe_metrics(yy[test],model.predict_proba(xx[test])[:,1]); rows.append({"seed":seed,"strategy":"property_matched_cached_decoys",**m}); props.append({"seed":seed,"selected_count":len(selected),**{f"{c}_mean":selected[c].mean() for c in COLUMNS}})
    metrics=pd.DataFrame(rows); metrics.to_csv(output/"decoy_subset_metrics.csv",index=False); pd.DataFrame(props).to_csv(output/"decoy_subset_properties.csv",index=False); manifest={"seeds":list(seeds),"cached_pool":"data/raw/xai/dude_hivpr/decoys_final.ism","new_downloads":False,"property_matching":"3000-seed candidate pool; closest 1000 to active standardized descriptor centroid","split":"scaffold-grouped 70/15/15","models_are_validation_only":True}; write_json(output/"decoy_subset_manifest.json",manifest); summary={"pr_auc_std":float(metrics.pr_auc.std(ddof=1)),"pr_auc_mean":float(metrics.pr_auc.mean()),"pr_auc_min":float(metrics.pr_auc.min()),"pr_auc_max":float(metrics.pr_auc.max()),"instability_warning":bool(metrics.pr_auc.std(ddof=1)>=.05),"limitation":"DUD-E decoys are computational and not experimentally verified inactive compounds."}; write_json(output/"decoy_sensitivity_summary.json",summary); return metrics
