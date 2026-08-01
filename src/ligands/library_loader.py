"""Load small screening libraries from CSV, SMI, or SDF inputs."""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Any

from rdkit import Chem

logger = logging.getLogger(__name__)


class LibraryLoaderError(RuntimeError):
    """Raised when a screening library cannot be loaded."""


def load_library(input_path: Path, max_compounds: int = 500, id_column: str = "compound_id", smiles_column: str = "smiles") -> dict[str, Any]:
    """Load a library from CSV, SMI or SDF and return an import summary."""

    if not input_path.exists():
        raise LibraryLoaderError(f"Library input not found: {input_path}")

    records: list[dict[str, Any]] = []
    if input_path.suffix.lower() == ".csv":
        with input_path.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        for row in rows:
            compound_id = row.get(id_column, "").strip()
            smiles = row.get(smiles_column, "").strip()
            if not compound_id:
                raise LibraryLoaderError("Every library row requires a compound_id")
            if not smiles:
                raise LibraryLoaderError(f"Compound {compound_id} is missing a SMILES value")
            records.append({"compound_id": compound_id, "smiles": smiles, "source": row.get("source", "csv")})
    elif input_path.suffix.lower() == ".smi":
        for line in input_path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) < 2:
                raise LibraryLoaderError(f"SMI line is malformed: {line}")
            compound_id, smiles = parts[0], parts[1]
            records.append({"compound_id": compound_id, "smiles": smiles, "source": "smi"})
    elif input_path.suffix.lower() == ".sdf":
        supplier = Chem.SDMolSupplier(str(input_path), removeHs=False, sanitize=True)
        for index, mol in enumerate(supplier, start=1):
            if mol is None:
                continue
            compound_id = mol.GetProp("compound_id") if mol.HasProp("compound_id") else f"mol_{index:03d}"
            smiles = Chem.MolToSmiles(mol)
            records.append({"compound_id": compound_id, "smiles": smiles, "source": "sdf"})
    else:
        raise LibraryLoaderError(f"Unsupported library format: {input_path.suffix}")

    if len(records) > max_compounds:
        raise LibraryLoaderError(f"Library size {len(records)} exceeds the configured maximum {max_compounds}")

    summary = {"input_path": str(input_path), "record_count": len(records), "records": records[:max_compounds]}
    return summary
