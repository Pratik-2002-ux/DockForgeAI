from __future__ import annotations
import numpy as np
def local_top_contributions(values,feature_names,count=10):
 order=np.argsort(-np.abs(values))[:count]; return [{"feature":feature_names[i],"contribution":float(values[i])} for i in order]
