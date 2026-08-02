"""Ensemble and conformal uncertainty without silent biological claims."""
from __future__ import annotations
import numpy as np
def ensemble_uncertainty(probabilities):
 values=np.asarray(probabilities,float); return {"mean":values.mean(0),"standard_deviation":values.std(0),"entropy":-(values.mean(0)*np.log(values.mean(0)+1e-12)+(1-values.mean(0))*np.log(1-values.mean(0)+1e-12))}
def conformal_threshold(y_cal,p_cal,alpha=.1):
 scores=np.where(np.asarray(y_cal)==1,1-np.asarray(p_cal),np.asarray(p_cal)); quantile=min(1.,np.ceil((len(scores)+1)*(1-alpha))/len(scores)); return float(np.quantile(scores,quantile,method="higher"))
def conformal_sets(p,threshold): return [(0 if value<=threshold else None,1 if 1-value<=threshold else None) for value in np.asarray(p)]
