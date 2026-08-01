"""Resolve structure duplicates and active/decoy overlap."""

from __future__ import annotations

from typing import Any


def resolve_duplicates(records: list[dict[str, Any]], active_keys: set[str] | None = None) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Deduplicate by canonical SMILES, full InChIKey, then audit connectivity."""
    active_keys = active_keys or set(); accepted, excluded = [], []; smiles_seen = set(); key_seen = set(); connectivity: dict[str, str] = {}
    for record in records:
        reason = None; smiles = record["canonical_isomeric_smiles"]; key = record["inchikey"]; block = record["inchikey_connectivity"]
        if key in active_keys: reason = "overlaps_active_structure"
        elif smiles in smiles_seen: reason = "duplicate_canonical_smiles"
        elif key in key_seen: reason = "duplicate_inchikey"
        if reason:
            excluded.append({**record, "exclusion_reason": reason}); continue
        if block in connectivity and connectivity[block] != key: record = {**record, "connectivity_duplicate_group": block}
        connectivity[block] = key; smiles_seen.add(smiles); key_seen.add(key); accepted.append(record)
    return accepted, excluded
