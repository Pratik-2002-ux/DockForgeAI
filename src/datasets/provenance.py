"""Dataset provenance and checksum utilities."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def sha256_file(path: Path) -> str:
    """Return a file SHA-256 digest."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""): digest.update(block)
    return digest.hexdigest()


def retrieval_timestamp() -> str:
    """Return an ISO-8601 UTC timestamp."""
    return datetime.now(timezone.utc).isoformat()


def write_provenance(path: Path, sources: list[dict[str, Any]], release: str | None = None, retrieved_at: str | None = None) -> Path:
    """Write source provenance without private absolute paths."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {"retrieved_at": retrieved_at or retrieval_timestamp(), "source_release": release or "not reported by source", "sources": sources}
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return path
