"""Select a diverse shortlist from the ranked hits."""

from __future__ import annotations

import logging
from typing import Any

logger = logging.getLogger(__name__)


def select_diverse_hits(rows: list[dict[str, Any]], top_n: int = 10, similarity_threshold: float = 0.7) -> list[dict[str, Any]]:
    """Greedily select a diverse set of hits while preserving the highest-ranked compound."""

    selected: list[dict[str, Any]] = []
    for row in sorted(rows, key=lambda item: item.get("rank", 0)):
        if not selected:
            selected.append(row)
            continue
        if all(abs(float(row.get("consensus_score", 0.0)) - float(item.get("consensus_score", 0.0))) > 0.0 for item in selected):
            selected.append(row)
        if len(selected) >= top_n:
            break
    return selected
