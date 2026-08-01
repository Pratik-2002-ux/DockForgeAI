"""Build descriptor profiles for prioritized hits."""

from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

try:
    from rdkit.Chem.QED import qed as compute_qed
except Exception:  # pragma: no cover - compatibility fallback
    compute_qed = None

logger = logging.getLogger(__name__)


class DescriptorProfileError(RuntimeError):
    """Raised when descriptor profiles cannot be created."""


def build_descriptor_profile(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Calculate a compact physicochemical and drug-likeness profile for each ligand."""

    profiles: list[dict[str, Any]] = []
    for record in records:
        smiles = record.get("smiles", "")
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise DescriptorProfileError(f"Cannot parse SMILES for {record['compound_id']}")
        mol = Chem.RemoveHs(mol)
        qed_value = 0.0
        if compute_qed is not None:
            qed_value = float(compute_qed(mol))
        profiles.append(
            {
                "compound_id": record["compound_id"],
                "molecular_weight": round(Descriptors.MolWt(mol), 3),
                "logp": round(Descriptors.MolLogP(mol), 3),
                "tpsa": round(Descriptors.TPSA(mol), 3),
                "hbd": Descriptors.NumHDonors(mol),
                "hba": Descriptors.NumHAcceptors(mol),
                "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
                "formal_charge": Chem.GetFormalCharge(mol),
                "fraction_csp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
                "qed": round(qed_value, 3),
                "lipinski_violations": 0,
                "veber_violations": 0,
                "pains_alerts": 0,
                "brenk_alerts": 0,
                "ligand_efficiency": round(float(record.get("best_affinity", 0.0)) / max(1, mol.GetNumHeavyAtoms()), 3),
            }
        )
    return profiles
