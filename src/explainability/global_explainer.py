from __future__ import annotations
import numpy as np
def model_importance(model,feature_names):
 estimator=model.steps[-1][1] if hasattr(model,"steps") else model
 if hasattr(estimator,"coef_"): values=np.asarray(estimator.coef_[0])
 elif hasattr(estimator,"feature_importances_"): values=np.asarray(estimator.feature_importances_)
 else: raise ValueError("Model has no native global importance")
 return sorted([{"feature":name,"importance":float(abs(value)),"direction":float(np.sign(value))} for name,value in zip(feature_names,values)],key=lambda row:row["importance"],reverse=True)
