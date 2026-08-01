import numpy as np
from src.ml.applicability_domain import applicability_domain
def test_domain_output_fields():
 result=applicability_domain(np.array([[1,0],[0,1]],bool),np.array([[1,0]],bool),np.array([[0.],[1.]]),np.array([[0.1]]),["a","b"])[0];assert "in_domain" in result and result["nearest_training_compound_id"]=="a"
def test_domain_confidence_bounded():
 assert 0<=applicability_domain(np.array([[1,0],[0,1]],bool),np.array([[1,0]],bool),np.array([[0.],[1.]]),np.array([[0.1]]),["a","b"])[0]["domain_confidence_score"]<=1
