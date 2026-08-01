"""Compute a transparent consensus ranking for prioritized hits."""

from __future__ import annotations

import logging
from typing import Any

logger = logging.getLogger(__name__)


class ConsensusRankerError(RuntimeError):
    """Raised when consensus ranking cannot be computed."""


def compute_consensus_ranking(rows: list[dict[str, Any]], weights: dict[str, float]) -> list[dict[str, Any]]:
    """Compute a weighted consensus score from docking, interactions, drug-likeness and efficiency features."""

    total_weight = sum(weights.values())
    if abs(total_weight - 1.0) > 1e-9:
        raise ConsensusRankerError("Weights must sum to 1.0")

    ranked: list[dict[str, Any]] = []
    for row in rows:
        score = 0.0
        score += weights.get("docking", 0.0) * float(row.get("norm_best_affinity", 0.0))
        score += weights.get("interaction", 0.0) * float(row.get("norm_interaction", 0.0))
        score += weights.get("druglikeness", 0.0) * float(row.get("norm_qed", 0.0))
        score += weights.get("structural_alerts", 0.0) * float(row.get("norm_alerts", 0.0))
        score += weights.get("ligand_efficiency", 0.0) * float(row.get("norm_ligand_efficiency", 0.0))
        ranked.append({**row, "consensus_score": round(score, 4)})
    ranked.sort(key=lambda item: item["consensus_score"], reverse=True)
    for index, item in enumerate(ranked, start=1):
        item["rank"] = index
    return ranked
