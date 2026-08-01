"""Sequential, cached ChEMBL activity API client."""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

import requests


class ChEMBLClientError(RuntimeError):
    """Raised when ChEMBL data cannot be retrieved or validated."""


class ChEMBLClient:
    """Retrieve target activities page-by-page with retry and resume."""

    endpoint = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
    status_endpoint = "https://www.ebi.ac.uk/chembl/api/data/status.json"

    def __init__(self, cache_dir: Path, timeout: int = 60, retries: int = 4, page_size: int = 1000, session: Any | None = None) -> None:
        self.cache_dir, self.timeout, self.retries, self.page_size = cache_dir, timeout, retries, page_size
        self.session = session or requests.Session()
        self.session.headers.update({"User-Agent": "DockForgeAI-XAI-Dataset/1.0 (public scientific dataset curation)"})

    def _request(self, params: dict[str, Any]) -> dict[str, Any]:
        last: Exception | None = None
        for attempt in range(self.retries + 1):
            try:
                response = self.session.get(self.endpoint, params=params, timeout=self.timeout)
                response.raise_for_status(); payload = response.json()
                if not isinstance(payload.get("activities"), list): raise ChEMBLClientError("ChEMBL response lacks activities")
                return payload
            except (requests.RequestException, ValueError, ChEMBLClientError) as exc:
                last = exc
                if attempt == self.retries: break
                time.sleep(min(2**attempt, 8))
        raise ChEMBLClientError(f"ChEMBL request failed: {last}")

    def fetch_activities(self, target_id: str = "CHEMBL243", refresh: bool = False, offline: bool = False) -> list[dict[str, Any]]:
        """Fetch all activities, caching each validated page."""
        self.cache_dir.mkdir(parents=True, exist_ok=True); records = []; offset = 0
        while True:
            cache = self.cache_dir / f"activity_page_{offset:07d}.json"
            if cache.is_file() and not refresh:
                try: payload = json.loads(cache.read_text(encoding="utf-8"))
                except (OSError, json.JSONDecodeError) as exc: raise ChEMBLClientError(f"Invalid cache {cache.name}: {exc}") from exc
            elif offline:
                if offset == 0: raise ChEMBLClientError("Offline mode requested but no ChEMBL cache is available")
                break
            else:
                payload = self._request({"target_chembl_id": target_id, "limit": self.page_size, "offset": offset})
                cache.write_text(json.dumps(payload), encoding="utf-8")
            page = payload.get("activities", [])
            if not isinstance(page, list): raise ChEMBLClientError(f"Invalid activities in {cache.name}")
            records.extend(page)
            meta = payload.get("page_meta", {})
            next_url = meta.get("next")
            if not next_url or not page: break
            offset += self.page_size
        return records

    def fetch_status(self, offline: bool = False) -> dict[str, Any]:
        """Return cached ChEMBL service status/release metadata when available."""
        cache = self.cache_dir / "chembl_status.json"
        if cache.is_file():
            try: return json.loads(cache.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError): return {}
        if offline: return {}
        try:
            response = self.session.get(self.status_endpoint, timeout=self.timeout); response.raise_for_status(); payload = response.json()
            cache.write_text(json.dumps(payload), encoding="utf-8"); return payload
        except (requests.RequestException, ValueError, OSError): return {}
