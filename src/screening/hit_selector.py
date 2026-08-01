"""Select refinement candidates from screening scores."""

from __future__ import annotations

import csv
import logging
import math
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


class HitSelectorError(RuntimeError):
    """Raised when refinement hits cannot be selected."""


def select_hits(score_table: Path, output_path: Path, top_fraction: float = 0.10, minimum_count: int = 10, maximum_count: int = 30) -> list[dict[str, Any]]:
    """Select the top fraction of successful ligands as refinement candidates."""

    with score_table.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))

    successful = [row for row in rows if row.get("status") == "success"]
    successful = sorted(successful, key=lambda item: float(item.get("best_affinity", 0.0)))
    count = max(minimum_count, math.ceil(len(successful) * top_fraction))
    count = min(count, maximum_count, len(successful))
    selected = successful[:count]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["compound_id", "best_affinity", "status", "runtime_seconds", "output_path"])
        writer.writeheader()
        writer.writerows(selected)
    return selected
