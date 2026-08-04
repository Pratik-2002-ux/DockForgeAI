"""Frozen Morgan-Tanimoto applicability-domain calculation."""
from __future__ import annotations

import numpy as np


def applicability_domain(query: np.ndarray, reference: np.ndarray, ids, pki, threshold: float) -> dict:
    q = query.astype(np.int32, copy=False)
    r = reference.astype(np.int32, copy=False)
    intersection = q @ r.T
    union = q.sum(1)[:, None] + r.sum(1)[None, :] - intersection
    similarity = np.divide(intersection, union, out=np.zeros_like(intersection, dtype=float), where=union > 0)
    nearest = similarity.argmax(axis=1)
    top_n = min(5, similarity.shape[1])
    top5 = np.partition(similarity, similarity.shape[1] - top_n, axis=1)[:, -top_n:].mean(axis=1)
    maximum = similarity[np.arange(len(similarity)), nearest]
    return {
        "maximum_training_tanimoto": maximum,
        "mean_top5_tanimoto": top5,
        "nearest_training_compound_id": np.asarray(ids)[nearest],
        "nearest_training_compound_pki": np.asarray(pki)[nearest],
        "descriptor_space_distance": np.full(len(query), np.nan),
        "applicability_domain_status": np.where(maximum >= threshold, "in_domain", "out_of_domain"),
    }
