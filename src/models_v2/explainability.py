"""Model-behaviour explanations for v2 predictions."""
from __future__ import annotations
import numpy as np
from sklearn.inspection import permutation_importance
def similarity_explanation(query_similarity,neighbor_ids,neighbor_labels): return {"maximum_similarity":float(np.max(query_similarity)),"neighbor_ids":list(neighbor_ids),"neighbor_labels":list(neighbor_labels),"statement":"Similarity is model context, not causal biological evidence."}
def permutation_explanation(model,x,y,n_repeats=5,seed=2026): return permutation_importance(model,x,y,n_repeats=n_repeats,random_state=seed,scoring="average_precision").importances_mean
def counterfactual_status(): return {"available":False,"reason":"Molecule-generating counterfactual optimization is intentionally not implemented in the core v2 predictor."}
