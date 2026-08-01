"""Evaluate redocking success against the crystallographic reference."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[2]


class RedockingValidationError(RuntimeError):
    """Raised when the redocking analysis cannot be summarized."""


def validate_redocking(
    analysis_csv: Path,
    summary_path: Path,
    threshold: float,
) -> dict[str, Any]:
    """Summarize redocking performance from the RMSD analysis CSV."""

    if not analysis_csv.exists():
        raise RedockingValidationError(f"Analysis CSV not found: {analysis_csv}")

    with analysis_csv.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise RedockingValidationError("No docking poses were analyzed")

    parsed_rows: list[dict[str, Any]] = []
    for row in rows:
        parsed_rows.append(
            {
                "pose_rank": int(row["pose_rank"]),
                "affinity_kcal_mol": float(row["affinity_kcal_mol"]),
                "vina_rmsd_lower_bound_angstrom": float(row["vina_rmsd_lower_bound_angstrom"]),
                "vina_rmsd_upper_bound_angstrom": float(row["vina_rmsd_upper_bound_angstrom"]),
                "crystal_in_place_rmsd_angstrom": float(row["crystal_in_place_rmsd_angstrom"]),
                "aligned_best_rmsd_angstrom": float(row["aligned_best_rmsd_angstrom"]),
                "passes_2A_threshold": row["passes_2A_threshold"].lower() == "true",
                "atom_records": int(row["atom_records"]),
                "pose_pdbqt": row["pose_pdbqt"],
            }
        )

    sorted_rows = sorted(parsed_rows, key=lambda item: item["pose_rank"])
    best_affinity_row = min(sorted_rows, key=lambda item: item["affinity_kcal_mol"])
    minimum_rmsd_row = min(sorted_rows, key=lambda item: item["crystal_in_place_rmsd_angstrom"])

    passing_rows = [row for row in sorted_rows if row["crystal_in_place_rmsd_angstrom"] <= threshold]
    first_passing_rank = next((row["pose_rank"] for row in sorted_rows if row["crystal_in_place_rmsd_angstrom"] <= threshold), None)

    overall_success = any(row["crystal_in_place_rmsd_angstrom"] <= threshold for row in sorted_rows)
    top_rank_success = sorted_rows[0]["crystal_in_place_rmsd_angstrom"] <= threshold

    summary = {
        "pose_count": len(sorted_rows),
        "best_affinity_kcal_mol": round(best_affinity_row["affinity_kcal_mol"], 3),
        "best_affinity_pose_rank": best_affinity_row["pose_rank"],
        "rmsd_best_affinity_pose_angstrom": round(best_affinity_row["crystal_in_place_rmsd_angstrom"], 3),
        "min_rmsd_angstrom": round(minimum_rmsd_row["crystal_in_place_rmsd_angstrom"], 3),
        "min_rmsd_pose_rank": minimum_rmsd_row["pose_rank"],
        "affinity_min_rmsd_pose_kcal_mol": round(minimum_rmsd_row["affinity_kcal_mol"], 3),
        "first_passing_pose_rank": first_passing_rank,
        "passing_pose_count": len(passing_rows),
        "passing_pose_percentage": round((len(passing_rows) / len(sorted_rows)) * 100.0, 1) if sorted_rows else 0.0,
        "overall_pose_reproduction_success": overall_success,
        "top_ranked_pose_success": top_rank_success,
        "threshold_angstrom": threshold,
        "interpretation": (
            "Validated strongly" if top_rank_success else "Pose reproduction observed but ranking imperfect" if overall_success else "Failed pose reproduction"
        ),
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary
