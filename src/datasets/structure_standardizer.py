"""RDKit structure standardization and descriptor calculation."""

from __future__ import annotations

import hashlib
import math
from typing import Any

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.Scaffolds import MurckoScaffold


class StructureStandardizationError(ValueError):
    """Raised when a chemical structure is invalid."""


def _catalog(name: str) -> FilterCatalog:
    params = FilterCatalogParams()
    params.AddCatalog(getattr(FilterCatalogParams.FilterCatalogs, name))
    return FilterCatalog(params)


_PAINS = _catalog("PAINS")
_BRENK = _catalog("BRENK")


def _largest_organic_fragment(molecule: Chem.Mol) -> Chem.Mol:
    fragments = [fragment for fragment in Chem.GetMolFrags(molecule, asMols=True, sanitizeFrags=True) if any(atom.GetAtomicNum() == 6 for atom in fragment.GetAtoms())]
    if not fragments: raise StructureStandardizationError("No organic fragment")
    return max(fragments, key=lambda mol: (mol.GetNumHeavyAtoms(), Descriptors.MolWt(mol)))


def standardize_structure(smiles: str) -> dict[str, Any]:
    """Standardize a SMILES structure while preserving stereochemistry."""
    if not isinstance(smiles, str) or not smiles.strip(): raise StructureStandardizationError("Missing SMILES")
    molecule = Chem.MolFromSmiles(smiles, sanitize=True)
    if molecule is None: raise StructureStandardizationError("Invalid SMILES")
    try:
        molecule = _largest_organic_fragment(molecule); Chem.SanitizeMol(molecule)
        canonical = Chem.MolToSmiles(molecule, isomericSmiles=True, canonical=True)
        molecule = Chem.MolFromSmiles(canonical)
        inchi = Chem.MolToInchi(molecule); inchikey = Chem.InchiToInchiKey(inchi)
    except Exception as exc:
        raise StructureStandardizationError(f"Sanitization failed: {exc}") from exc
    if not inchikey: raise StructureStandardizationError("InChIKey generation failed")
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=molecule, includeChirality=True) or canonical
    values = {
        "canonical_isomeric_smiles": canonical, "inchi": inchi, "inchikey": inchikey,
        "inchikey_connectivity": inchikey.split("-")[0],
        "molecular_formula": rdMolDescriptors.CalcMolFormula(molecule),
        "molecular_weight": Descriptors.MolWt(molecule), "exact_molecular_weight": Descriptors.ExactMolWt(molecule),
        "logp": Crippen.MolLogP(molecule), "tpsa": rdMolDescriptors.CalcTPSA(molecule),
        "hbd": Lipinski.NumHDonors(molecule), "hba": Lipinski.NumHAcceptors(molecule),
        "rotatable_bonds": Lipinski.NumRotatableBonds(molecule), "ring_count": rdMolDescriptors.CalcNumRings(molecule),
        "aromatic_ring_count": rdMolDescriptors.CalcNumAromaticRings(molecule), "formal_charge": Chem.GetFormalCharge(molecule),
        "heavy_atom_count": molecule.GetNumHeavyAtoms(), "fraction_csp3": rdMolDescriptors.CalcFractionCSP3(molecule),
        "qed": QED.qed(molecule), "pains_alert_count": len(_PAINS.GetMatches(molecule)),
        "brenk_alert_count": len(_BRENK.GetMatches(molecule)), "murcko_scaffold_smiles": scaffold,
        "scaffold_id": "SCF-" + hashlib.sha256(scaffold.encode()).hexdigest()[:16],
    }
    for key, value in values.items():
        if isinstance(value, float) and not math.isfinite(value): raise StructureStandardizationError(f"Non-finite descriptor: {key}")
    return values
