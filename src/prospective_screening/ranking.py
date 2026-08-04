"""Deterministic, review-only ranking without docking information."""
from __future__ import annotations

import numpy as np
import pandas as pd


DOMAIN_ORDER = {"in_domain": 1, "borderline": 2, "out_of_domain": 3}


def rank_predictions(frame: pd.DataFrame) -> pd.DataFrame:
    ranked = frame.copy()
    ranked["predicted_pki_rank"] = ranked["predicted_pki"].rank(method="first", ascending=False).astype(int)
    ranked["uncertainty_rank"] = pd.Series(pd.NA, index=ranked.index, dtype="Int64")
    if ranked["prediction_uncertainty"].notna().any():
        ranked["uncertainty_rank"] = ranked["prediction_uncertainty"].rank(method="first", ascending=True).astype("Int64")
    ranked["domain_rank"] = ranked["applicability_domain_status"].map(DOMAIN_ORDER).astype(int)
    ordering = ranked.sort_values(["domain_rank", "predicted_pki", "compound_id"], ascending=[True, False, True])
    scaffold_first = ordering.drop_duplicates("murcko_scaffold_smiles").index
    ranked["diversity_rank"] = pd.Series(np.arange(1, len(scaffold_first) + 1), index=scaffold_first).reindex(ranked.index).astype("Int64")
    novelty = 1.0 - ranked["maximum_training_tanimoto"]
    alerts = ranked[["pains_flag", "brenk_flag"]].astype(int).sum(axis=1)
    domain_score = ranked["applicability_domain_status"].map({"in_domain": 2.0, "borderline": 1.0, "out_of_domain": 0.0})
    pki_scaled = ranked["predicted_pki"].rank(pct=True)
    ranked["review_priority_score"] = domain_score + pki_scaled + 0.25 * novelty - 0.15 * alerts
    ranked["review_priority_rank"] = ranked.sort_values(["review_priority_score", "predicted_pki", "compound_id"], ascending=[False, False, True]).index.to_series().rank(method="first").reindex(ranked.index).astype(int)
    # Replace the index-based intermediate rank with an exact sequential order.
    order = ranked.sort_values(["review_priority_score", "predicted_pki", "compound_id"], ascending=[False, False, True]).index
    ranked.loc[order, "review_priority_rank"] = np.arange(1, len(ranked) + 1)
    return ranked


def review_pools(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    ordered = frame.sort_values("review_priority_rank")
    top100 = ordered.head(100).copy()
    diverse = ordered.drop_duplicates("murcko_scaffold_smiles").head(30).copy()
    if len(diverse) < min(30, len(ordered)):
        diverse = pd.concat([diverse, ordered.loc[~ordered.index.isin(diverse.index)]])[:30]
    return top100, diverse
