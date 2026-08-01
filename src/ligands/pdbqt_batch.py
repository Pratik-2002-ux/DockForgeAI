"""Prepare ligand PDBQT files in a small batch workflow."""

from __future__ import annotations

import json
import logging
import shutil
import subprocess
from pathlib import Path
from typing import Any

import yaml

logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SETTINGS_PATH = PROJECT_ROOT / "config" / "settings.yaml"


class PDBQTBatchError(RuntimeError):
    """Raised when ligand PDBQT generation fails."""


def load_settings(settings_path: Path = SETTINGS_PATH) -> dict[str, Any]:
    """Load project settings from YAML."""

    with settings_path.open("r", encoding="utf-8") as handle:
        settings = yaml.safe_load(handle)
    if not isinstance(settings, dict):
        raise ValueError("Settings file must contain a YAML dictionary")
    return settings


def prepare_pdbqt(records: list[dict[str, Any]], output_dir: Path, overwrite: bool = False) -> dict[str, Any]:
    """Generate ligand PDBQT files using Meeko when available."""

    settings = load_settings()
    output_dir.mkdir(parents=True, exist_ok=True)
    meeko_executable = settings.get("docking", {}).get("meeko_export_executable", "mk_export.py")
    meeko_path = shutil.which(meeko_executable) or meeko_executable

    prepared: list[dict[str, Any]] = []
    failed: list[dict[str, Any]] = []
    for record in records:
        sdf_path = Path(record["output_path"])
        pdbqt_path = output_dir / f"{record['compound_id']}.pdbqt"
        if pdbqt_path.exists() and not overwrite:
            prepared.append({"compound_id": record["compound_id"], "output_path": str(pdbqt_path), "status": "skipped"})
            continue
        command = [meeko_path, "-s", str(pdbqt_path), str(sdf_path)]
        completed = subprocess.run(command, cwd=str(PROJECT_ROOT), capture_output=True, text=True, check=False)
        if completed.returncode == 0 and pdbqt_path.exists() and pdbqt_path.stat().st_size > 0:
            text = pdbqt_path.read_text(encoding="utf-8")
            if "ROOT" in text and "ENDROOT" in text and "TORSDOF" in text:
                prepared.append({"compound_id": record["compound_id"], "output_path": str(pdbqt_path), "status": "ready"})
            else:
                failed.append({"compound_id": record["compound_id"], "reason": "missing_pdbqt_torsion_markers"})
        else:
            failed.append({"compound_id": record["compound_id"], "reason": "meeko_failed", "stderr": completed.stderr})
    return {"prepared": prepared, "failed": failed}
