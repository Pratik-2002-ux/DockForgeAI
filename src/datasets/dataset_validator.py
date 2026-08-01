"""Hard validation rules for the final HIVPR XAI dataset."""

from __future__ import annotations

import hashlib
import math
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem

from src.datasets.scaffold_splitter import scaffold_leakage

MANDATORY = ("compound_id", "standardized_compound_id", "source_database", "source_record_id", "source_molecule_id", "source_target_id", "original_smiles", "canonical_isomeric_smiles", "inchikey", "binary_label", "label_source", "murcko_scaffold_smiles", "scaffold_id", "dataset_split", "retrieval_date", "source_release", "provenance_reference")
FINITE = ("molecular_weight", "exact_molecular_weight", "logp", "tpsa", "hbd", "hba", "rotatable_bonds", "formal_charge", "heavy_atom_count", "qed")


class DatasetValidationError(RuntimeError):
    """Raised when one or more hard dataset rules fail."""


def validate_records(records: list[dict[str, Any]]) -> dict[str, Any]:
    """Validate identifiers, structures, labels, splits, provenance, and leakage."""
    errors = []
    for column in MANDATORY:
        if any(row.get(column) in (None, "") for row in records): errors.append(f"missing mandatory field: {column}")
    smiles = [row.get("canonical_isomeric_smiles") for row in records]
    if len(smiles) != len(set(smiles)): errors.append("duplicate standardized structures")
    active = {row["inchikey"] for row in records if row.get("binary_label") == "active"}; decoy = {row["inchikey"] for row in records if row.get("binary_label") == "inactive_or_decoy"}
    if active & decoy: errors.append("active/decoy structure overlap")
    if any(row.get("binary_label") not in {"active", "inactive_or_decoy"} for row in records): errors.append("invalid binary label")
    if any(row.get("dataset_split") not in {"train", "validation", "test"} for row in records): errors.append("invalid dataset split")
    leakage = scaffold_leakage(records)
    if leakage: errors.append(f"scaffold leakage: {len(leakage)} scaffolds")
    if any(Path(str(value)).is_absolute() for row in records for value in row.values() if isinstance(value, str)): errors.append("absolute path in dataset")
    for row in records:
        for column in FINITE:
            try:
                if not math.isfinite(float(row[column])): raise ValueError
            except (KeyError, TypeError, ValueError): errors.append(f"non-finite {column}: {row.get('compound_id')}")
    if errors: raise DatasetValidationError("; ".join(sorted(set(errors))))
    return {"row_count": len(records), "active_count": len(active), "decoy_count": len(decoy), "scaffold_count": len({row["scaffold_id"] for row in records}), "scaffold_leakage_count": 0}


def validate_output_files(csv_path: Path, parquet_path: Path, sdf_path: Path) -> dict[str, int]:
    """Ensure CSV, Parquet and SDF contain identical row counts."""
    csv_count = len(pd.read_csv(csv_path)); parquet_count = len(pd.read_parquet(parquet_path))
    sdf_count = sum(molecule is not None for molecule in Chem.SDMolSupplier(str(sdf_path), removeHs=False))
    if len({csv_count, parquet_count, sdf_count}) != 1: raise DatasetValidationError(f"Output row-count mismatch: CSV={csv_count}, Parquet={parquet_count}, SDF={sdf_count}")
    return {"csv_count": csv_count, "parquet_count": parquet_count, "sdf_count": sdf_count}
