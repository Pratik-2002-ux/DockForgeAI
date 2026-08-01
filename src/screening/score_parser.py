"""Parse screening docking outputs into a score table."""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


class ScoreParserError(RuntimeError):
    """Raised when docking outputs cannot be parsed."""


def parse_scores(results: list[dict[str, Any]], output_path: Path) -> list[dict[str, Any]]:
    """Convert batch results into a simple CSV table."""

    rows: list[dict[str, Any]] = []
    for result in results:
        rows.append(
            {
                "compound_id": result["compound_id"],
                "best_affinity": result.get("best_affinity", None),
                "mean_top3_affinity": result.get("mean_top3_affinity", None),
                "median_affinity": result.get("median_affinity", None),
                "pose_count": result.get("pose_count", 0),
                "status": result.get("status", "failed"),
                "runtime_seconds": result.get("runtime_seconds", None),
                "error_message": result.get("error_message", ""),
                "output_path": result.get("output_path", ""),
            }
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["compound_id", "best_affinity", "mean_top3_affinity", "median_affinity", "pose_count", "status", "runtime_seconds", "error_message", "output_path"])
        writer.writeheader()
        writer.writerows(rows)
    return rows
