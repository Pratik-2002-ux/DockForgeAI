"""Official DUD-E HIVPR structure-file client."""

from __future__ import annotations

import hashlib
import time
from pathlib import Path
from typing import Any

import requests


class DUDEClientError(RuntimeError):
    """Raised when DUD-E files cannot be retrieved or parsed."""


class DUDEClient:
    """Download and inspect official DUD-E target files."""

    base_url = "http://dude.docking.org/targets"
    filenames = ("decoys_final.ism", "actives_final.ism", "receptor.pdb", "crystal_ligand.mol2")

    def __init__(self, cache_dir: Path, timeout: int = 60, retries: int = 4, session: Any | None = None) -> None:
        self.cache_dir, self.timeout, self.retries = cache_dir, timeout, retries
        self.session = session or requests.Session(); self.session.headers.update({"User-Agent": "DockForgeAI-XAI-Dataset/1.0"})

    def download(self, target_id: str = "HIVPR", refresh: bool = False, offline: bool = False) -> dict[str, Path]:
        """Download known files individually; the decoy file is mandatory."""
        self.cache_dir.mkdir(parents=True, exist_ok=True); outputs = {}
        for name in self.filenames:
            path = self.cache_dir / name
            if path.is_file() and path.stat().st_size and not refresh: outputs[name] = path; continue
            if offline: continue
            url = f"{self.base_url}/{target_id.lower()}/{name}"; last = None
            for attempt in range(self.retries + 1):
                try:
                    response = self.session.get(url, timeout=self.timeout); response.raise_for_status()
                    if not response.content: raise DUDEClientError(f"Empty DUD-E response for {name}")
                    path.write_bytes(response.content); outputs[name] = path; break
                except (requests.RequestException, OSError, DUDEClientError) as exc:
                    last = exc
                    if attempt < self.retries: time.sleep(min(2**attempt, 8))
            if name == "decoys_final.ism" and name not in outputs: raise DUDEClientError(f"DUD-E decoys unavailable: {last}")
        if "decoys_final.ism" not in outputs: raise DUDEClientError("Offline mode requested but decoys_final.ism is unavailable")
        return outputs

    @staticmethod
    def parse_ism(path: Path) -> list[dict[str, str]]:
        """Parse an inspected DUD-E whitespace-delimited SMILES/ID file."""
        records = []
        for line_number, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), start=1):
            if not line.strip() or line.startswith("#"): continue
            parts = line.split()
            if len(parts) < 2: raise DUDEClientError(f"Malformed {path.name} line {line_number}")
            records.append({"original_smiles": parts[0], "source_molecule_id": parts[1], "source_record_id": f"{path.name}:{line_number}"})
        return records

    @staticmethod
    def checksum(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()
