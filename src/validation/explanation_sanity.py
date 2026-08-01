"""Non-causal explanation sanity controls."""
from __future__ import annotations
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from .core import OUTPUT, load_features, safe_metrics, write_json

def run_explanation_sanity(output=OUTPUT/"explanation_sanity"):
    x,y,index=load_features("descriptor"); tr=index.dataset_split.eq("train").to_numpy(); te=index.dataset_split.eq("test").to_numpy(); model=LogisticRegression(max_iter=1000,random_state=2026).fit(x[tr],y[tr]); base=safe_metrics(y[te],model.predict_proba(x[te])[:,1]); native=np.abs(model.coef_[0]); perm=permutation_importance(model,x[te],y[te],n_repeats=5,random_state=2026,scoring="average_precision").importances_mean; agreement=float(spearmanr(native,perm).statistic)
    order=np.argsort(native); removal=[]
    for k in (10,20,30):
        for kind,drop in (("top",order[-k:]),("bottom",order[:k])):
            keep=np.setdiff1d(np.arange(x.shape[1]),drop); candidate=LogisticRegression(max_iter=1000,random_state=2026).fit(x[tr][:,keep],y[tr]); m=safe_metrics(y[te],candidate.predict_proba(x[te][:,keep])[:,1]); removal.append({"k":k,"removed":kind,"pr_auc":m["pr_auc"],"pr_auc_loss":base["pr_auc"]-m["pr_auc"]})
    rng=np.random.default_rng(2026); randomized=rng.permutation(model.coef_.ravel()); random_corr=float(spearmanr(model.coef_.ravel(),randomized).statistic); labels=rng.permutation(y[tr]); scrambled=LogisticRegression(max_iter=1000,random_state=2026).fit(x[tr],labels); label_corr=float(spearmanr(np.abs(model.coef_[0]),np.abs(scrambled.coef_[0])).statistic)
    pd.DataFrame([{"test":"parameter_randomization","importance_correlation":random_corr},{"test":"label_randomization","importance_correlation":label_corr}]).to_csv(output/"randomization_tests.csv",index=False); pd.DataFrame(removal).to_csv(output/"feature_removal_tests.csv",index=False); pd.DataFrame([{"comparison":"coefficient_vs_permutation_importance","spearman_correlation":agreement}]).to_csv(output/"importance_agreement.csv",index=False)
    logits=model.decision_function(x[te]); reconstructed=x[te]@model.coef_.ravel()+model.intercept_[0]; fidelity=float(np.max(np.abs(logits-reconstructed))); pd.DataFrame([{"model":"descriptor_logistic_regression","maximum_reconstruction_error":fidelity}]).to_csv(output/"local_fidelity.csv",index=False)
    summary={"explanations_are_causal":False,"statement":"Explanations describe model behaviour and do not establish causal chemical mechanisms.","importance_agreement":agreement,"parameter_randomization_correlation":random_corr,"label_randomization_correlation":label_corr,"maximum_local_fidelity_error":fidelity,"feature_removal":removal}; write_json(output/"explanation_sanity_summary.json",summary); return summary
