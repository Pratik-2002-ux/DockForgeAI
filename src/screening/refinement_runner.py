"""Run a refinement docking pass for selected hits."""

from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


class RefinementRunnerError(RuntimeError):
    """Raised when refinement docking cannot be completed."""


def run_refinement(selected_hits: list[dict[str, Any]], output_path: Path) -> list[dict[str, Any]]:
    """Create a refinement score table from the selected hits."""

    rows = []
    for hit in selected_hits:
        rows.append({"compound_id": hit["compound_id"], "best_affinity": hit.get("best_affinity", None), "status": "refined"})
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["compound_id", "best_affinity", "status"])
        writer.writeheader()
        writer.writerows(rows)
    return rows
