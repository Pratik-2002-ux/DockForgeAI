import numpy as np
from src.explainability.local_explainer import local_top_contributions
def test_local_contributions_rank_absolute(): assert local_top_contributions(np.array([1.,-3.]),["a","b"])[0]["feature"]=="b"
def test_local_contribution_count(): assert len(local_top_contributions(np.array([1.,2.]),["a","b"],1))==1
