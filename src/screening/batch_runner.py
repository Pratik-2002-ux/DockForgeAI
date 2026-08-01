"""Run docking for a batch of ligands sequentially with resume support."""

from __future__ import annotations

import csv
import json
import logging
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any

import yaml

logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SETTINGS_PATH = PROJECT_ROOT / "config" / "settings.yaml"


class BatchRunnerError(RuntimeError):
    """Raised when a docking batch cannot be executed."""


def load_settings(settings_path: Path = SETTINGS_PATH) -> dict[str, Any]:
    """Load project settings from YAML."""

    with settings_path.open("r", encoding="utf-8") as handle:
        settings = yaml.safe_load(handle)
    if not isinstance(settings, dict):
        raise ValueError("Settings file must contain a YAML dictionary")
    return settings


def run_batch(records: list[dict[str, Any]], config_template: Path, output_dir: Path, timeout_seconds: int = 1800, resume: bool = True) -> dict[str, Any]:
    """Run docking for each ligand using Vina and capture summary data."""

    settings = load_settings()
    output_dir.mkdir(parents=True, exist_ok=True)
    vina = settings.get("docking", {}).get("vina_executable", "vina")
    vina_path = shutil.which(vina) or vina
    results: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    for record in records:
        compound_id = record["compound_id"]
        output_path = output_dir / compound_id / "output.pdbqt"
        if output_path.exists() and resume:
            results.append({"compound_id": compound_id, "status": "skipped", "output_path": str(output_path)})
            continue
        config_path = output_dir / compound_id / "config.txt"
        config_path.parent.mkdir(parents=True, exist_ok=True)
        config_text = config_template.read_text(encoding="utf-8").replace("{ligand}", str(PROJECT_ROOT / "data/processed/screening_library/pdbqt" / f"{compound_id}.pdbqt"))
        config_path.write_text(config_text, encoding="utf-8")
        start = time.perf_counter()
        completed = subprocess.run([vina_path, "--config", str(config_path)], capture_output=True, text=True, check=False, timeout=timeout_seconds)
        runtime = round(time.perf_counter() - start, 3)
        if completed.returncode == 0 and output_path.exists():
            results.append({"compound_id": compound_id, "status": "success", "runtime_seconds": runtime, "output_path": str(output_path)})
        else:
            failures.append({"compound_id": compound_id, "status": "failed", "runtime_seconds": runtime, "stderr": completed.stderr})
    return {"results": results, "failures": failures}
