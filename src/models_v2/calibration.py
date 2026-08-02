"""Calibration metrics and temperature scaling."""
from __future__ import annotations
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import expit,logit
from sklearn.metrics import brier_score_loss
def ece(y,p,bins=10):
 y=np.asarray(y); p=np.asarray(p); edges=np.linspace(0,1,bins+1); return float(sum(mask.mean()*abs(y[mask].mean()-p[mask].mean()) for lo,hi in zip(edges[:-1],edges[1:]) if (mask:=((p>=lo)&(p<(hi if hi<1 else hi+1e-12)))).any()))
def temperature_scale(y,p):
 y=np.asarray(y); logits=logit(np.clip(p,1e-7,1-1e-7)); loss=lambda t:float(-np.mean(y*np.log(expit(logits/t)+1e-12)+(1-y)*np.log(1-expit(logits/t)+1e-12))); temperature=float(minimize_scalar(loss,bounds=(.05,10),method="bounded").x); return temperature,expit(logits/temperature)
def calibration_metrics(y,p): return {"ece":ece(y,p),"brier_score":float(brier_score_loss(y,p))}
