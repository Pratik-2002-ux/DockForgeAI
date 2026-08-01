"""Safe report generation and download helpers for Streamlit."""

from __future__ import annotations

import io
import json
import subprocess
import sys
import threading
import zipfile
from pathlib import Path
from typing import Any, Iterable

PROJECT_ROOT = Path(__file__).resolve().parents[2]
_GENERATION_LOCK = threading.Lock()


def format_file_size(size_bytes: int) -> str:
    """Format a non-negative byte count for display."""

    if size_bytes < 0: raise ValueError("File size cannot be negative")
    size = float(size_bytes)
    for unit in ("B", "KB", "MB", "GB"):
        if size < 1024 or unit == "GB": return f"{size:.0f} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} GB"


def safe_download_files(candidates: Iterable[Path], project_root: Path = PROJECT_ROOT) -> list[Path]:
    """Return existing regular files inside the project, without duplicates."""

    root = project_root.resolve(); accepted: list[Path] = []; seen: set[Path] = set()
    for candidate in candidates:
        path = candidate.resolve()
        try: path.relative_to(root)
        except ValueError: continue
        if path.is_file() and path not in seen: accepted.append(path); seen.add(path)
    return accepted


def build_zip_bytes(files: Iterable[Path], project_root: Path = PROJECT_ROOT) -> bytes:
    """Build an in-memory ZIP using project-relative archive names."""

    root = project_root.resolve(); buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in safe_download_files(files, root): archive.write(path, path.relative_to(root).as_posix())
    return buffer.getvalue()


def load_manifest(project_root: Path = PROJECT_ROOT) -> dict[str, Any]:
    """Load the report manifest without exposing host paths."""

    path = project_root / "data/processed/final_reports/report_manifest.json"
    try:
        value = json.loads(path.read_text(encoding="utf-8")); return value if isinstance(value, dict) else {}
    except (OSError, json.JSONDecodeError): return {}


def generate_reports(timeout_seconds: int = 180, overwrite: bool = True, project_root: Path = PROJECT_ROOT) -> tuple[bool, str]:
    """Run only the report CLI with a timeout and duplicate-request lock."""

    if not _GENERATION_LOCK.acquire(blocking=False): return False, "Report generation is already in progress."
    try:
        command = [sys.executable, str(project_root / "scripts/generate_final_report.py")]
        if overwrite: command.append("--overwrite")
        result = subprocess.run(command, cwd=project_root, capture_output=True, text=True, timeout=timeout_seconds, check=False)
        message = (result.stdout if result.returncode == 0 else result.stderr).strip()
        return result.returncode == 0, message or "Report command completed."
    except subprocess.TimeoutExpired: return False, f"Report generation exceeded the {timeout_seconds}-second timeout."
    except OSError as exc: return False, f"Report generation could not start: {exc}"
    finally: _GENERATION_LOCK.release()


def categorized_exports(project_root: Path = PROJECT_ROOT) -> dict[str, list[Path]]:
    """Collect fixed, safe export categories without accepting user paths."""

    root = project_root
    categories = {
        "Reports": list((root / "data/processed/final_reports").glob("**/*")),
        "Tables": [root / path for path in (
            "data/processed/docking_tables/1HSG_MK1_redocking_analysis.csv", "data/processed/interaction_reports/interactions.csv",
            "data/processed/screening_tables/screening_scores.csv", "data/processed/prioritization_reports/consensus_ranking.csv")],
        "Summaries": list((root / "data/processed").glob("**/*summary.json")) + [root / "data/processed/screening_reports/screening_report.json", root / "data/processed/prioritization_reports/hit_prioritization_report.json"],
        "Molecule images": list((root / "data/processed/report_assets/molecules").glob("*")),
        "Figures": list((root / "data/processed/report_assets/figures").glob("*")),
    }
    return {name: safe_download_files(files, root) for name, files in categories.items()}
