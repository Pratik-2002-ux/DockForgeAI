from __future__ import annotations
import numpy as np
from scipy.stats import pearsonr,spearmanr,kendalltau
from sklearn.metrics import mean_squared_error,mean_absolute_error,median_absolute_error,r2_score
def metrics(y,p):
 y=np.asarray(y);p=np.asarray(p);err=np.abs(y-p)
 def corr(fn):
  try:return float(fn(y,p).statistic)
  except:return None
 return {"rmse":float(mean_squared_error(y,p)**.5),"mae":float(mean_absolute_error(y,p)),"median_ae":float(median_absolute_error(y,p)),"r2":float(r2_score(y,p)),"pearson":corr(pearsonr),"spearman":corr(spearmanr),"kendall":corr(kendalltau),"within_0_5":float((err<=.5).mean()),"within_1_0":float((err<=1).mean())}
