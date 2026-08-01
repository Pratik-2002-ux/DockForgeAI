"""Generate hit-prioritization reports and tables."""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


class PrioritizationReportError(RuntimeError):
    """Raised when prioritization reports cannot be generated."""


def generate_prioritization_report(rows: list[dict[str, Any]], output_dir: Path) -> dict[str, Any]:
    """Write ranking tables and reports for prioritized hits."""

    output_dir.mkdir(parents=True, exist_ok=True)
    ranking_csv = output_dir / "consensus_ranking.csv"
    diverse_csv = output_dir / "diverse_top_hits.csv"
    md_path = output_dir / "hit_prioritization_report.md"
    json_path = output_dir / "hit_prioritization_report.json"

    normalized_rows = []
    for row in rows:
        normalized = {key: row.get(key) for key in ["compound_id", "rank", "consensus_score", "best_affinity", "qed", "pains_alerts", "brenk_alerts", "alerts", "interaction", "ligand_efficiency", "smiles"]}
        normalized_rows.append(normalized)

    with ranking_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["compound_id", "rank", "consensus_score", "best_affinity", "qed", "pains_alerts", "brenk_alerts", "alerts", "interaction", "ligand_efficiency", "smiles"])
        writer.writeheader()
        writer.writerows(normalized_rows)
    with diverse_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["compound_id", "rank", "consensus_score", "best_affinity", "qed", "pains_alerts", "brenk_alerts", "alerts", "interaction", "ligand_efficiency", "smiles"])
        writer.writeheader()
        writer.writerows([{key: row.get(key) for key in ["compound_id", "rank", "consensus_score", "best_affinity", "qed", "pains_alerts", "brenk_alerts", "alerts", "interaction", "ligand_efficiency", "smiles"]} for row in rows[:5]])
    md_path.write_text("# Hit prioritization report\n\n" + json.dumps(rows[:5], indent=2) + "\n", encoding="utf-8")
    json_path.write_text(json.dumps({"rows": rows}, indent=2), encoding="utf-8")
    return {"ranking_csv": str(ranking_csv), "diverse_csv": str(diverse_csv), "md_path": str(md_path), "json_path": str(json_path)}
