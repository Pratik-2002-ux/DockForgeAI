"""Deterministic no-replacement descriptor matching."""
from __future__ import annotations
import numpy as np, pandas as pd
from sklearn.preprocessing import StandardScaler
from .core import OUTPUT, best_predictions, load_data, safe_metrics, write_json

MATCH_COLUMNS=["molecular_weight","logp","tpsa","hbd","hba","rotatable_bonds","formal_charge","heavy_atom_count"]
def standardized_mean_difference(a,b):
    den=np.sqrt((np.var(a,ddof=1)+np.var(b,ddof=1))/2); return float((np.mean(a)-np.mean(b))/den) if den else 0.
def run_descriptor_matching(caliper=.2, output=OUTPUT/"matched_subsets"):
    d=load_data(); test=d[d.dataset_split.eq("test")].copy(); active=test[test.binary_label.eq("active")]; decoy=test[~test.binary_label.eq("active")]; scaler=StandardScaler().fit(test[MATCH_COLUMNS]); xa=scaler.transform(active[MATCH_COLUMNS]); xd=scaler.transform(decoy[MATCH_COLUMNS]); available=set(range(len(decoy))); pairs=[]
    for i,row in enumerate(xa):
        if not available: break
        candidates=np.array(sorted(available)); distances=np.linalg.norm(xd[candidates]-row,axis=1)/np.sqrt(len(MATCH_COLUMNS)); j=int(candidates[np.argmin(distances)]); distance=float(np.min(distances))
        if distance<=caliper: pairs.append({"active_id":active.iloc[i].compound_id,"decoy_id":decoy.iloc[j].compound_id,"distance":distance}); available.remove(j)
    pair=pd.DataFrame(pairs,columns=["active_id","decoy_id","distance"]); pair.to_csv(output/"matched_pairs.csv",index=False); matched=set(pair.active_id)|set(pair.decoy_id); test[~test.compound_id.isin(matched)][["compound_id","binary_label"]].to_csv(output/"unmatched_compounds.csv",index=False)
    def balance(frame):
        aa=frame[frame.binary_label.eq("active")]; bb=frame[~frame.binary_label.eq("active")]; return pd.DataFrame([{"descriptor":c,"standardized_mean_difference":standardized_mean_difference(aa[c],bb[c]),"absolute_smd":abs(standardized_mean_difference(aa[c],bb[c]))} for c in MATCH_COLUMNS])
    before=balance(test); after=balance(test[test.compound_id.isin(matched)]) if matched else before.copy(); before.to_csv(output/"matching_balance_before.csv",index=False); after.to_csv(output/"matching_balance_after.csv",index=False)
    pred=best_predictions().set_index("compound_id"); full=safe_metrics((test.binary_label=="active").astype(int),test.compound_id.map(pred.predicted_probability)); subset=test[test.compound_id.isin(matched)]; matched_metrics=safe_metrics((subset.binary_label=="active").astype(int),subset.compound_id.map(pred.predicted_probability)) if len(subset) else {}
    summary={"caliper":caliper,"replacement":False,"matched_pairs":len(pair),"mean_absolute_smd_before":float(before.absolute_smd.mean()),"mean_absolute_smd_after":float(after.absolute_smd.mean()) if matched else None,"full_test":full,"matched_test":matched_metrics,"pr_auc_drop":float(full["pr_auc"]-matched_metrics["pr_auc"]) if matched_metrics else None,"status":"evaluated" if matched_metrics else "not_evaluated_no_pairs_within_caliper"}; write_json(output/"matched_subset_metrics.json",summary); return summary
