"""Generate interaction analysis reports."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


class InteractionReportError(RuntimeError):
    """Raised when interaction reports cannot be generated."""


def generate_interaction_report(rows: list[dict[str, Any]], output_dir: Path) -> dict[str, Any]:
    """Write Markdown, JSON and CSV outputs for interaction analysis."""

    output_dir.mkdir(parents=True, exist_ok=True)
    report_md_path = output_dir / "report.md"
    report_json_path = output_dir / "report.json"
    interactions_csv_path = output_dir / "interactions.csv"

    with interactions_csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["interaction_type", "protein_chain", "protein_residue_name", "protein_residue_number", "protein_atom", "ligand_atom", "distance_angstrom", "angle_degrees", "source_engine", "pose_label", "pose_rank", "notes"])
        writer.writeheader()
        writer.writerows(rows)

    report_md_path.write_text("# Interaction report\n\n" + "\n".join([f"- {row['interaction_type']} {row['protein_residue_name']}{row['protein_residue_number']}" for row in rows]) + "\n", encoding="utf-8")
    report_json_path.write_text(json.dumps({"row_count": len(rows), "rows": rows}, indent=2), encoding="utf-8")
    return {"report_md": str(report_md_path), "report_json": str(report_json_path), "interactions_csv": str(interactions_csv_path)}
