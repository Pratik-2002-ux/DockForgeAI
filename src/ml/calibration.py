from __future__ import annotations
import numpy as np
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression
from src.ml.metrics import classification_metrics
class ProbabilityCalibrator:
 def __init__(self,method="none"): self.method=method; self.model=None
 def fit(self,p,y):
  if self.method=="sigmoid": self.model=LogisticRegression().fit(np.asarray(p).reshape(-1,1),y)
  elif self.method=="isotonic": self.model=IsotonicRegression(out_of_bounds="clip").fit(p,y)
  return self
 def predict(self,p):
  if self.method=="none": return np.asarray(p)
  if self.method=="sigmoid": return self.model.predict_proba(np.asarray(p).reshape(-1,1))[:,1]
  return self.model.predict(p)
def choose_calibration(validation_y,validation_probabilities):
 candidates=[]
 for method in ("none","sigmoid","isotonic"):
  cal=ProbabilityCalibrator(method).fit(validation_probabilities,validation_y); probs=cal.predict(validation_probabilities); metrics=classification_metrics(validation_y,probs); candidates.append((metrics["brier_score"]+metrics["expected_calibration_error"],method,cal,metrics))
 return min(candidates,key=lambda x:x[0])[1:]
