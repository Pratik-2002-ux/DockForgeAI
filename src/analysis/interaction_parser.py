"""Standardize interaction records into a common schema."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


class InteractionParserError(RuntimeError):
    """Raised when interaction data cannot be parsed."""


def parse_interactions(records: list[dict[str, Any]], output_dir: Path, pose_label: str, pose_rank: int) -> dict[str, Any]:
    """Write standardized interaction records to CSV/JSON."""

    output_dir.mkdir(parents=True, exist_ok=True)
    standardized: list[dict[str, Any]] = []
    for record in records:
        standardized.append(
            {
                "interaction_type": record.get("interaction_type", "unknown"),
                "protein_chain": record.get("protein_chain", ""),
                "protein_residue_name": record.get("protein_residue_name", ""),
                "protein_residue_number": int(record.get("protein_residue_number", 0)),
                "protein_atom": record.get("protein_atom", ""),
                "ligand_atom": record.get("ligand_atom", ""),
                "distance_angstrom": record.get("distance_angstrom"),
                "angle_degrees": record.get("angle_degrees"),
                "source_engine": record.get("source_engine", "fallback"),
                "pose_label": pose_label,
                "pose_rank": pose_rank,
                "notes": record.get("notes", ""),
            }
        )

    csv_path = output_dir / f"{pose_label}_interactions.csv"
    json_path = output_dir / f"{pose_label}_interactions.json"

    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "interaction_type",
                "protein_chain",
                "protein_residue_name",
                "protein_residue_number",
                "protein_atom",
                "ligand_atom",
                "distance_angstrom",
                "angle_degrees",
                "source_engine",
                "pose_label",
                "pose_rank",
                "notes",
            ],
        )
        writer.writeheader()
        writer.writerows(standardized)

    json_path.write_text(json.dumps(standardized, indent=2), encoding="utf-8")
    return {"pose_label": pose_label, "pose_rank": pose_rank, "row_count": len(standardized), "csv_path": str(csv_path), "json_path": str(json_path)}
