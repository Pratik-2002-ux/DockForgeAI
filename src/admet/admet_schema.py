"""Validate optional external ADMET input schemas."""

from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


class ADMETSchemaError(RuntimeError):
    """Raised when an ADMET input table does not conform to the expected schema."""


def validate_admet_table(path: Path) -> dict[str, Any]:
    """Validate that the external ADMET table contains the required columns and valid values."""

    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise ADMETSchemaError("External ADMET table is empty")
    if "compound_id" not in rows[0] or "source_tool" not in rows[0]:
        raise ADMETSchemaError("External ADMET table must include compound_id and source_tool")

    for row in rows:
        compound_id = row.get("compound_id", "").strip()
        if not compound_id:
            raise ADMETSchemaError("Every ADMET row needs a compound_id")
        source_tool = row.get("source_tool", "").strip()
        if not source_tool:
            raise ADMETSchemaError(f"Compound {compound_id} is missing source_tool")
        for key in ["ames_probability", "herg_probability", "hepatotoxicity_probability", "bbb_probability", "gi_absorption_probability", "cyp3a4_inhibition_probability", "cyp2d6_inhibition_probability"]:
            value = row.get(key, "")
            if value in {"", None}:
                continue
            try:
                score = float(value)
            except ValueError as exc:
                raise ADMETSchemaError(f"Invalid numeric value for {key} in {compound_id}") from exc
            if not 0.0 <= score <= 1.0:
                raise ADMETSchemaError(f"Probability for {compound_id} is outside [0, 1]")
    return {"row_count": len(rows)}
