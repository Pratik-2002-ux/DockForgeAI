from __future__ import annotations
import numpy as np
def linear_contributions(model,x,feature_names):
 estimator=model.steps[-1][1] if hasattr(model,"steps") else model; values=x.toarray() if hasattr(x,"toarray") else np.asarray(x)
 coefficients=np.asarray(estimator.coef_[0]); contributions=values*coefficients
 return contributions,feature_names
