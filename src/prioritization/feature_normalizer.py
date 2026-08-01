"""Normalize features used for consensus ranking."""

from __future__ import annotations

import logging
from typing import Any

logger = logging.getLogger(__name__)


def normalize_features(rows: list[dict[str, Any]], beneficial: set[str], risk: set[str]) -> list[dict[str, Any]]:
    """Normalize features using a simple min–max transform."""

    normalized: list[dict[str, Any]] = []
    for feature in beneficial | risk:
        values = [float(row.get(feature, 0.0)) for row in rows if row.get(feature) is not None]
        if not values:
            continue
        min_value = min(values)
        max_value = max(values)
        if max_value == min_value:
            for row in rows:
                row[f"norm_{feature}"] = 1.0
            continue
        for row in rows:
            value = float(row.get(feature, 0.0))
            if feature in beneficial:
                row[f"norm_{feature}"] = (value - min_value) / (max_value - min_value)
            else:
                row[f"norm_{feature}"] = (max_value - value) / (max_value - min_value)
    return rows
