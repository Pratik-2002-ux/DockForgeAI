"""Generate reports for the virtual-screening workflow."""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


class ScreeningReportError(RuntimeError):
    """Raised when screening reports cannot be generated."""


def generate_screening_report(summary: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    """Write Markdown and JSON reports for the screening workflow."""

    output_dir.mkdir(parents=True, exist_ok=True)
    md_path = output_dir / "screening_report.md"
    json_path = output_dir / "screening_report.json"
    md_path.write_text("# Screening report\n\n" + json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return {"md_path": str(md_path), "json_path": str(json_path)}
