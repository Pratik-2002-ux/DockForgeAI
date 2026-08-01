"""Compare interaction patterns across poses."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


class InteractionComparatorError(RuntimeError):
    """Raised when interaction comparison fails."""


def _residue_key(record: dict[str, Any]) -> str:
    chain = record.get("protein_chain", "") or "A"
    residue_name = record.get("protein_residue_name", "") or ""
    residue_number = record.get("protein_residue_number", 0) or 0
    return f"{chain}:{residue_number}"


def compare_interactions(crystal: list[dict[str, Any]], docked: list[dict[str, Any]], output_path: Path) -> dict[str, Any]:
    """Compare crystallographic and docked interaction residues using a simple fingerprint."""

    crystal_residues = {_residue_key(record) for record in crystal}
    docked_residues = {_residue_key(record) for record in docked}

    shared = sorted(crystal_residues)
    lost = sorted(crystal_residues - docked_residues)
    gained = sorted(docked_residues - crystal_residues)

    union = crystal_residues | docked_residues
    jaccard_similarity = len(shared) / len(union) if union else 1.0

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps({"shared_residue_contacts": shared, "lost_residue_contacts": lost, "gained_residue_contacts": gained, "jaccard_similarity": jaccard_similarity}, indent=2), encoding="utf-8")
    return {
        "shared_residue_contacts": shared,
        "lost_residue_contacts": lost,
        "gained_residue_contacts": gained,
        "jaccard_similarity": jaccard_similarity,
    }
