"""Parse Vina multi-model PDBQT output into per-pose files and scoring rows."""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[2]


class ResultParsingError(RuntimeError):
    """Raised when Vina output cannot be parsed."""


def parse_vina_result(
    vina_output_path: Path,
    output_dir: Path,
    score_table_path: Path,
    metadata_path: Path,
) -> dict[str, Any]:
    """Parse a multi-model Vina PDBQT file into separate pose files and CSV rows."""

    if not vina_output_path.exists():
        raise ResultParsingError(f"Vina output not found: {vina_output_path}")

    output_dir.mkdir(parents=True, exist_ok=True)
    score_table_path.parent.mkdir(parents=True, exist_ok=True)

    text = vina_output_path.read_text(encoding="utf-8")
    model_blocks = re.findall(r"MODEL\s+(\d+)(.*?)(?:ENDMDL\s*)", text, flags=re.DOTALL)
    if not model_blocks:
        raise ResultParsingError("No MODEL blocks found in the Vina output")

    pose_rows: list[dict[str, Any]] = []
    pose_paths: list[Path] = []

    for index, (model_id, block) in enumerate(model_blocks, start=1):
        result_line = next(
            (line for line in block.splitlines() if line.startswith("REMARK VINA RESULT")),
            None,
        )
        if result_line is None:
            raise ResultParsingError(f"Pose {index} is missing a VINA RESULT remark")

        parts = result_line.split()
        if len(parts) < 6:
            raise ResultParsingError(f"Pose {index} has an unexpected VINA RESULT line")
        affinity = float(parts[-3])
        rmsd_lower = float(parts[-2])
        rmsd_upper = float(parts[-1])

        atom_records = sum(1 for line in block.splitlines() if line.startswith(("ATOM", "HETATM")))

        pose_path = output_dir / f"pose_{index:03d}.pdbqt"
        pose_path.write_text(f"MODEL {model_id}\n{block}ENDMDL\n", encoding="utf-8")
        pose_paths.append(pose_path)

        pose_rows.append(
            {
                "pose_rank": index,
                "affinity_kcal_mol": affinity,
                "vina_rmsd_lower_bound_angstrom": rmsd_lower,
                "vina_rmsd_upper_bound_angstrom": rmsd_upper,
                "atom_records": atom_records,
                "pose_pdbqt": str(pose_path),
            }
        )

    if [row["pose_rank"] for row in pose_rows] != list(range(1, len(pose_rows) + 1)):
        raise ResultParsingError("Parsed pose ranks are not sequential")

    with score_table_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "pose_rank",
                "affinity_kcal_mol",
                "vina_rmsd_lower_bound_angstrom",
                "vina_rmsd_upper_bound_angstrom",
                "atom_records",
                "pose_pdbqt",
            ],
        )
        writer.writeheader()
        writer.writerows(pose_rows)

    metadata = {
        "input_pdbqt": str(vina_output_path),
        "poses_written": len(pose_rows),
        "score_table": str(score_table_path),
        "pose_paths": [str(path) for path in pose_paths],
    }
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return metadata
