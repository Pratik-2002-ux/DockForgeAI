"""Ligand chemistry validation utilities for DockForgeAI."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import yaml
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SETTINGS_PATH = PROJECT_ROOT / "config" / "settings.yaml"


class LigandValidationError(RuntimeError):
    """Raised when a ligand cannot be read or validated."""


def load_settings(
    settings_path: Path = SETTINGS_PATH,
) -> dict[str, Any]:
    """Load project settings from YAML."""

    if not settings_path.exists():
        raise FileNotFoundError(
            f"Settings file not found: {settings_path}"
        )

    with settings_path.open("r", encoding="utf-8") as file:
        settings = yaml.safe_load(file)

    if not isinstance(settings, dict):
        raise ValueError(
            "Settings file must contain a YAML dictionary."
        )

    return settings


def load_sdf_molecule(sdf_path: Path) -> Chem.Mol:
    """Load the first valid molecule from an SDF file."""

    if not sdf_path.exists():
        raise FileNotFoundError(
            f"Ligand SDF file not found: {sdf_path}"
        )

    supplier = Chem.SDMolSupplier(
        str(sdf_path),
        removeHs=False,
        sanitize=True,
    )

    molecule = next(
        (
            molecule
            for molecule in supplier
            if molecule is not None
        ),
        None,
    )

    if molecule is None:
        raise LigandValidationError(
            f"RDKit could not parse ligand: {sdf_path}"
        )

    return molecule


def calculate_ligand_properties(
    molecule: Chem.Mol,
) -> dict[str, Any]:
    """Calculate molecular and drug-like properties."""

    molecule_without_hydrogens = Chem.RemoveHs(molecule)

    return {
        "molecular_formula": (
            rdMolDescriptors.CalcMolFormula(molecule)
        ),
        "molecular_weight": round(
            Descriptors.MolWt(molecule),
            3,
        ),
        "exact_molecular_weight": round(
            Descriptors.ExactMolWt(molecule),
            3,
        ),
        "heavy_atoms": molecule.GetNumHeavyAtoms(),
        "total_atoms": molecule.GetNumAtoms(),
        "formal_charge": Chem.GetFormalCharge(molecule),
        "hydrogen_bond_donors": (
            Lipinski.NumHDonors(molecule)
        ),
        "hydrogen_bond_acceptors": (
            Lipinski.NumHAcceptors(molecule)
        ),
        "rotatable_bonds": (
            Lipinski.NumRotatableBonds(molecule)
        ),
        "rings": Lipinski.RingCount(molecule),
        "aromatic_rings": (
            rdMolDescriptors.CalcNumAromaticRings(molecule)
        ),
        "tpsa": round(
            rdMolDescriptors.CalcTPSA(molecule),
            3,
        ),
        "logp": round(
            Descriptors.MolLogP(molecule),
            3,
        ),
        "canonical_smiles": Chem.MolToSmiles(
            molecule_without_hydrogens,
            canonical=True,
        ),
    }


def validate_ligand(
    sdf_path: Path,
) -> dict[str, Any]:
    """Validate ligand structure and return a report."""

    molecule = load_sdf_molecule(sdf_path)

    try:
        Chem.SanitizeMol(molecule)
        sanitization_status = "passed"
    except Exception as error:
        sanitization_status = f"failed: {error}"

    conformer_count = molecule.GetNumConformers()
    has_3d_coordinates = False

    if conformer_count > 0:
        conformer = molecule.GetConformer()
        has_3d_coordinates = bool(conformer.Is3D())

    return {
        "source_file": str(sdf_path),
        "ligand_name": sdf_path.stem,
        "sanitization_status": sanitization_status,
        "number_of_conformers": conformer_count,
        "has_3d_coordinates": has_3d_coordinates,
        "properties": calculate_ligand_properties(molecule),
    }


def save_validation_report(
    report: dict[str, Any],
) -> Path:
    """Save the validation report as JSON."""

    settings = load_settings()

    output_directory = (
        PROJECT_ROOT
        / settings["directories"]["ligand_reports"]
    )

    output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    ligand_name = report["ligand_name"]
    output_path = (
        output_directory
        / f"{ligand_name}_validation.json"
    )

    with output_path.open("w", encoding="utf-8") as file:
        json.dump(
            report,
            file,
            indent=2,
        )

    return output_path