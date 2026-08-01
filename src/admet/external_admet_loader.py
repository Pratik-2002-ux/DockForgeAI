"""Load optional external ADMET predictions supplied by the user."""

from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


class ExternalADMETLoaderError(RuntimeError):
    """Raised when external ADMET data cannot be loaded."""


def load_external_admet(path: Path) -> list[dict[str, Any]]:
    """Load an external ADMET CSV table and return rows."""

    if not path.exists():
        raise ExternalADMETLoaderError(f"ADMET input file not found: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    return rows
