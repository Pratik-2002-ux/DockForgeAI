"""Locked external identity exclusion."""
from __future__ import annotations
from pathlib import Path
import pandas as pd

def locked_identities(path: str | Path) -> dict[str,set[str]]:
    frame=pd.read_csv(path); return {c:set(frame[c].dropna().astype(str)) if c in frame else set() for c in ("canonical_isomeric_smiles","inchikey","inchikey_connectivity","murcko_scaffold_smiles")}

def mark_external_overlap(frame: pd.DataFrame, locked: dict[str,set[str]]) -> pd.DataFrame:
    out=frame.copy(); out["external_exact_smiles_overlap"]=out.canonical_isomeric_smiles.isin(locked["canonical_isomeric_smiles"]); out["external_inchikey_overlap"]=out.inchikey.isin(locked["inchikey"]); out["external_connectivity_overlap"]=out.inchikey_connectivity.isin(locked["inchikey_connectivity"]); out["locked_external_scaffold_overlap"]=out.murcko_scaffold_smiles.isin(locked["murcko_scaffold_smiles"]); out["external_overlap_flag"]=out[["external_exact_smiles_overlap","external_inchikey_overlap","external_connectivity_overlap"]].any(axis=1); return out
