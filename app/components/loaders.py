"""Data-loading utilities for the DockForgeAI Streamlit application."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def resolve_path(relative_path: str | Path) -> Path:
    """Resolve a project-relative path."""

    return PROJECT_ROOT / Path(relative_path)


def load_json(relative_path: str | Path) -> dict[str, Any]:
    """Load a JSON file safely."""

    path = resolve_path(relative_path)

    if not path.exists():
        return {}

    try:
        with path.open("r", encoding="utf-8") as handle:
            data = json.load(handle)
    except (json.JSONDecodeError, OSError):
        return {}

    return data if isinstance(data, dict) else {}


def load_csv(relative_path: str | Path) -> pd.DataFrame:
    """Load a CSV file safely."""

    path = resolve_path(relative_path)

    if not path.exists():
        return pd.DataFrame()

    try:
        return pd.read_csv(path)
    except (pd.errors.ParserError, OSError, UnicodeDecodeError):
        return pd.DataFrame()


def load_text(relative_path: str | Path) -> str:
    """Load a text or Markdown file safely."""

    path = resolve_path(relative_path)

    if not path.exists():
        return ""

    try:
        return path.read_text(encoding="utf-8")
    except OSError:
        return ""


def file_status(relative_path: str | Path) -> dict[str, Any]:
    """Return basic information for a project file."""

    path = resolve_path(relative_path)

    return {
        "path": str(path),
        "exists": path.exists(),
        "size_bytes": path.stat().st_size if path.exists() else 0,
    }