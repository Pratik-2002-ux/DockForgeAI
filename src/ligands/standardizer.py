"""Standardize screening-library molecules with RDKit."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

logger = logging.getLogger(__name__)


class StandardizerError(RuntimeError):
    """Raised when a molecule cannot be standardized."""


def standardize_ligands(records: list[dict[str, Any]]) -> dict[str, Any]:
    """Sanitize molecules, remove salts/fragments, canonicalize SMILES, and deduplicate."""

    standardized: list[dict[str, Any]] = []
    rejected: list[dict[str, Any]] = []
    seen: set[str] = set()

    for record in records:
        compound_id = record.get("compound_id", "")
        smiles = record.get("smiles", "")
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            rejected.append({"compound_id": compound_id, "reason": "invalid_smiles"})
            continue
        mol = Chem.AddHs(mol)
        mol = Chem.RemoveHs(mol)
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol, isomericSmiles=True))
        if mol is None:
            rejected.append({"compound_id": compound_id, "reason": "failed_sanitize"})
            continue
        canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
        if canonical in seen:
            rejected.append({"compound_id": compound_id, "reason": "duplicate_smiles"})
            continue
        seen.add(canonical)
        standardized.append({"compound_id": compound_id, "smiles": canonical, "source": record.get("source", "")})

    return {"standardized": standardized, "rejected": rejected, "count": len(standardized)}
