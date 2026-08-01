"""Run AutoDock Vina and capture execution metadata."""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SETTINGS_PATH = PROJECT_ROOT / "config" / "settings.yaml"


class VinaExecutionError(RuntimeError):
    """Raised when Vina execution fails."""


def load_settings(settings_path: Path = SETTINGS_PATH) -> dict[str, Any]:
    """Load project settings from YAML."""

    if not settings_path.exists():
        raise FileNotFoundError(f"Settings file not found: {settings_path}")

    with settings_path.open("r", encoding="utf-8") as handle:
        settings = yaml.safe_load(handle)

    if not isinstance(settings, dict):
        raise ValueError("Settings file must contain a YAML dictionary.")

    return settings


def _resolve_executable(executable_name: str) -> str:
    """Resolve a command-line executable from the environment."""

    executable_path = shutil.which(executable_name)
    if executable_path is None:
        raise VinaExecutionError(f"Executable not found: {executable_name}")
    return executable_path


def _supports_log_flag(executable_path: str) -> bool:
    """Check whether the installed Vina supports the --log flag."""

    completed = subprocess.run(
        [executable_path, "--help"],
        capture_output=True,
        text=True,
        check=False,
        cwd=str(PROJECT_ROOT),
    )
    combined_output = (completed.stdout or "") + (completed.stderr or "")
    return "--log" in combined_output


def run_vina(
    config_path: Path,
    output_path: Path,
    log_dir: Path,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Run Vina using the configured docking file and capture execution metadata."""

    settings = load_settings()
    docking_settings = settings.get("docking", {})
    timeout_seconds = int(docking_settings.get("timeout_seconds", 7200))
    executable_name = str(docking_settings.get("vina_executable", "vina"))
    executable_path = _resolve_executable(executable_name)

    if output_path.exists() and output_path.stat().st_size > 0 and not overwrite:
        raise VinaExecutionError(f"Target output already exists and overwrite=False: {output_path}")

    log_dir.mkdir(parents=True, exist_ok=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    stdout_path = log_dir / "vina_stdout.txt"
    stderr_path = log_dir / "vina_stderr.txt"
    metadata_path = log_dir / "vina_execution.json"
    temp_output = output_path.with_suffix(output_path.suffix + ".tmp")

    if temp_output.exists():
        temp_output.unlink()
    if output_path.exists() and overwrite:
        output_path.unlink()

    version_output = subprocess.run(
        [executable_path, "--version"],
        capture_output=True,
        text=True,
        check=False,
        cwd=str(PROJECT_ROOT),
    )
    vina_version = (version_output.stdout or version_output.stderr or "unknown").strip()

    command = [executable_path, "--config", str(config_path), "--out", str(temp_output)]
    use_log_flag = _supports_log_flag(executable_path)
    if use_log_flag:
        command.extend(["--log", str(stdout_path)])

    env = os.environ.copy()
    env.update(
        {
            "OMP_NUM_THREADS": "4",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "VECLIB_MAXIMUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
        }
    )

    start_time = datetime.now(timezone.utc).isoformat()
    completed = subprocess.run(
        command,
        cwd=str(PROJECT_ROOT),
        capture_output=True,
        text=True,
        timeout=timeout_seconds,
        check=False,
        env=env,
    )
    end_time = datetime.now(timezone.utc).isoformat()

    stdout_path.write_text(completed.stdout or "", encoding="utf-8")
    stderr_path.write_text(completed.stderr or "", encoding="utf-8")

    if not use_log_flag:
        stdout_path.write_text(completed.stdout or "", encoding="utf-8")

    if completed.returncode != 0:
        raise VinaExecutionError(
            f"Vina exited with code {completed.returncode}: {completed.stderr.strip() or completed.stdout.strip()}"
        )

    if not temp_output.exists() or temp_output.stat().st_size == 0:
        raise VinaExecutionError("Vina did not produce an output PDBQT")

    output_text = temp_output.read_text(encoding="utf-8")
    if "MODEL" not in output_text:
        raise VinaExecutionError("Vina output does not contain any MODEL records")
    if "REMARK VINA RESULT" not in output_text:
        raise VinaExecutionError("Vina output does not contain VINA RESULT remarks")

    temp_output.replace(output_path)

    metadata = {
        "command": " ".join(command),
        "executable": executable_path,
        "vina_version": vina_version,
        "start_time_utc": start_time,
        "end_time_utc": end_time,
        "runtime_seconds": round(
            (datetime.fromisoformat(end_time) - datetime.fromisoformat(start_time)).total_seconds(),
            3,
        ),
        "return_code": completed.returncode,
        "stdout_path": str(stdout_path),
        "stderr_path": str(stderr_path),
        "output_path": str(output_path),
        "used_log_flag": use_log_flag,
    }
    metadata_path.write_text(__import__("json").dumps(metadata, indent=2), encoding="utf-8")
    return metadata
