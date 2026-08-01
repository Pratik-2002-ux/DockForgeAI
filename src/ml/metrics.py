from __future__ import annotations
import numpy as np
from sklearn.metrics import *
def expected_calibration_error(y,p,bins=10):
 edges=np.linspace(0,1,bins+1); total=len(y); value=0.
 for lo,hi in zip(edges[:-1],edges[1:]):
  mask=(p>=lo)&(p<(hi if hi<1 else hi+1e-9))
  if mask.any(): value+=mask.mean()*abs(y[mask].mean()-p[mask].mean())
 return float(value)
def classification_metrics(y,p):
 pred=(p>=.5).astype(int); tn,fp,fn,tp=confusion_matrix(y,pred,labels=[0,1]).ravel()
 return {"roc_auc":roc_auc_score(y,p),"pr_auc":average_precision_score(y,p),"accuracy":accuracy_score(y,pred),"balanced_accuracy":balanced_accuracy_score(y,pred),"precision":precision_score(y,pred,zero_division=0),"recall":recall_score(y,pred),"specificity":tn/(tn+fp) if tn+fp else 0,"sensitivity":tp/(tp+fn) if tp+fn else 0,"negative_predictive_value":tn/(tn+fn) if tn+fn else 0,"f1":f1_score(y,pred),"mcc":matthews_corrcoef(y,pred),"cohen_kappa":cohen_kappa_score(y,pred),"brier_score":brier_score_loss(y,p),"log_loss":log_loss(y,p,labels=[0,1]),"expected_calibration_error":expected_calibration_error(y,p),"tn":int(tn),"fp":int(fp),"fn":int(fn),"tp":int(tp)}
