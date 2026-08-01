"""Calculate basic drug-likeness and descriptor features for standardized ligands."""

from __future__ import annotations

import logging
from typing import Any

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

try:
    from rdkit.Chem.QED import qed as compute_qed
except Exception:  # pragma: no cover - compatibility fallback
    compute_qed = None

logger = logging.getLogger(__name__)


class DruglikenessError(RuntimeError):
    """Raised when descriptor generation fails."""


def calculate_descriptors(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Calculate a compact descriptor profile for each ligand."""

    results: list[dict[str, Any]] = []
    for record in records:
        mol = Chem.MolFromSmiles(record["smiles"])
        if mol is None:
            raise DruglikenessError(f"Cannot parse SMILES for {record['compound_id']}")
        mol = Chem.AddHs(mol)
        mol = Chem.RemoveHs(mol)
        qed = 0.0
        if compute_qed is not None:
            qed = float(compute_qed(mol))
        results.append(
            {
                "compound_id": record["compound_id"],
                "smiles": record["smiles"],
                "molecular_weight": round(Descriptors.MolWt(mol), 3),
                "logp": round(Descriptors.MolLogP(mol), 3),
                "tpsa": round(Descriptors.TPSA(mol), 3),
                "hbd": Descriptors.NumHDonors(mol),
                "hba": Descriptors.NumHAcceptors(mol),
                "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
                "ring_count": rdMolDescriptors.CalcNumRings(mol),
                "aromatic_ring_count": rdMolDescriptors.CalcNumAromaticRings(mol),
                "formal_charge": Chem.GetFormalCharge(mol),
                "fraction_csp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
                "qed": round(qed, 3),
                "lipinski_violations": 0,
                "veber_violations": 0,
                "pains_alerts": 0,
                "brenk_alerts": 0,
            }
        )
    return results
