"""Exact reconstruction of the frozen Track 1 feature representation."""
from __future__ import annotations

import hashlib
import json

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem


FEATURE_SCHEMA = {
    "feature_representation": "morgan_r3",
    "generator": "RDKit MorganGenerator",
    "radius": 3,
    "n_bits": 2048,
    "use_chirality": False,
    "feature_order": "bit indices 0..2047",
    "preprocessing": "none",
}


def schema_checksum() -> str:
    payload = json.dumps(FEATURE_SCHEMA, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def rebuild_morgan_r3(smiles) -> np.ndarray:
    generator = AllChem.GetMorganGenerator(radius=3, fpSize=2048)
    matrix = np.zeros((len(smiles), 2048), dtype=np.uint8)
    for row, value in enumerate(smiles):
        mol = Chem.MolFromSmiles(str(value))
        if mol is None:
            raise ValueError(f"Invalid standardized structure at row {row}")
        DataStructs.ConvertToNumpyArray(generator.GetFingerprint(mol), matrix[row])
    return matrix


def validate_feature_schema(model, representation: str) -> dict:
    matches = representation == "morgan_r3" and getattr(model, "n_features_in_", None) == 2048
    return {**FEATURE_SCHEMA, "model_n_features": getattr(model, "n_features_in_", None), "schema_checksum": schema_checksum(), "matches_frozen_model": matches}
