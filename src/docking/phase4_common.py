"""Shared deterministic provenance, checksum, and report utilities for Phase 4."""
from __future__ import annotations

import hashlib
import html
import json
import os
import platform
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of *path*."""

    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_json(value: Any, path: Path) -> None:
    """Write JSON atomically."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", dir=path.parent, delete=False, encoding="utf-8") as handle:
        temporary = Path(handle.name)
        json.dump(value, handle, indent=2, default=str)
        handle.write("\n")
    os.replace(temporary, path)


def provenance(root: Path, seed: int = 2026) -> dict[str, Any]:
    """Collect the reproducibility fields required by Phase 4."""

    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root, text=True, capture_output=True, check=True).stdout.strip()
    return {
        "generated_timestamp": datetime.now(timezone.utc).isoformat(),
        "git_commit": commit,
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "random_seed": seed,
    }


def write_reports(title: str, sections: dict[str, Any], html_path: Path, pdf_path: Path) -> None:
    """Create compact HTML and PDF reports from computed values."""

    html_path.parent.mkdir(parents=True, exist_ok=True)
    rows = "".join(f"<tr><th>{html.escape(str(key))}</th><td>{html.escape(str(value))}</td></tr>" for key, value in sections.items())
    document = f"<!doctype html><html><head><meta charset='utf-8'><title>{html.escape(title)}</title><style>body{{font-family:Arial;margin:2rem;color:#172033}}table{{border-collapse:collapse;width:100%}}th,td{{border:1px solid #ccd3df;padding:.55rem;text-align:left}}th{{background:#eef2f7}}</style></head><body><h1>{html.escape(title)}</h1><table>{rows}</table><p>Research workflow only; no experimental claims.</p></body></html>"
    html_path.write_text(document, encoding="utf-8")
    pdf = canvas.Canvas(str(pdf_path), pagesize=letter)
    width, height = letter
    pdf.setTitle(title)
    pdf.setFont("Helvetica-Bold", 16)
    pdf.drawString(54, height - 54, title)
    y = height - 84
    pdf.setFont("Helvetica", 9)
    for key, value in sections.items():
        line = f"{key}: {value}"
        for start in range(0, len(line), 105):
            if y < 54:
                pdf.showPage(); pdf.setFont("Helvetica", 9); y = height - 54
            pdf.drawString(54, y, line[start:start + 105]); y -= 13
    pdf.drawString(54, max(y - 12, 42), "Research workflow only; no experimental claims.")
    pdf.save()


def write_checksums(paths: list[Path], destination: Path) -> dict[str, str]:
    """Checksum existing outputs and write their mapping atomically."""

    checksums = {path.name: sha256(path) for path in paths}
    atomic_json(checksums, destination)
    return checksums
