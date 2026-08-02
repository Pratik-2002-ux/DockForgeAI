def replicate_cv(values):
 import numpy as np
 a=np.asarray(values,float);return float(a.std(ddof=1)/abs(a.mean())) if len(a)>1 and a.mean()!=0 else None
