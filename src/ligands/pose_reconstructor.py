"""Crystallographic ligand-pose reconstruction for DockForgeAI.

This module transfers chemically correct bond orders from an authoritative
SDF template onto ligand coordinates extracted from a protein PDB file.

The output preserves the crystallographic heavy-atom coordinates and is
suitable as the reference pose for later redocking RMSD calculations.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import yaml
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SETTINGS_PATH = PROJECT_ROOT / "config" / "settings.yaml"


class PoseReconstructionError(RuntimeError):
    """Raised when crystallographic pose reconstruction fails."""


def load_settings(
    settings_path: Path = SETTINGS_PATH,
) -> dict[str, Any]:
    """Load DockForgeAI settings."""

    if not settings_path.exists():
        raise FileNotFoundError(
            f"Settings file not found: {settings_path}"
        )

    with settings_path.open("r", encoding="utf-8") as handle:
        settings = yaml.safe_load(handle)

    if not isinstance(settings, dict):
        raise ValueError(
            "Settings file must contain a YAML dictionary."
        )

    return settings


def load_template_molecule(template_path: Path) -> Chem.Mol:
    """Load and sanitize a chemically correct SDF template."""

    if not template_path.exists():
        raise FileNotFoundError(
            f"Template SDF was not found: {template_path}"
        )

    supplier = Chem.SDMolSupplier(
        str(template_path),
        removeHs=False,
        sanitize=True,
    )

    template = next(
        (
            molecule
            for molecule in supplier
            if molecule is not None
        ),
        None,
    )

    if template is None:
        raise PoseReconstructionError(
            f"RDKit could not read template: {template_path}"
        )

    # The extracted crystallographic PDB contains heavy atoms only.
    template = Chem.RemoveHs(template)

    Chem.SanitizeMol(template)

    return template


def load_crystallographic_pose(pdb_path: Path) -> Chem.Mol:
    """Load an extracted PDB ligand without initially sanitizing it."""

    if not pdb_path.exists():
        raise FileNotFoundError(
            f"Extracted ligand PDB was not found: {pdb_path}"
        )

    pose = Chem.MolFromPDBFile(
        str(pdb_path),
        removeHs=False,
        sanitize=False,
        proximityBonding=True,
    )

    if pose is None:
        raise PoseReconstructionError(
            f"RDKit could not read crystallographic ligand: {pdb_path}"
        )

    if pose.GetNumAtoms() == 0:
        raise PoseReconstructionError(
            "The crystallographic ligand contains no atoms."
        )

    if pose.GetNumConformers() == 0:
        raise PoseReconstructionError(
            "The crystallographic ligand has no coordinates."
        )

    return pose


def get_coordinates(molecule: Chem.Mol) -> list[tuple[float, float, float]]:
    """Return coordinates from the first conformer."""

    if molecule.GetNumConformers() == 0:
        raise PoseReconstructionError(
            "Molecule contains no conformer."
        )

    conformer = molecule.GetConformer()

    return [
        (
            float(conformer.GetAtomPosition(index).x),
            float(conformer.GetAtomPosition(index).y),
            float(conformer.GetAtomPosition(index).z),
        )
        for index in range(molecule.GetNumAtoms())
    ]


def calculate_coordinate_difference(
    before: list[tuple[float, float, float]],
    after: list[tuple[float, float, float]],
) -> dict[str, float]:
    """Measure coordinate changes after bond-order assignment."""

    if len(before) != len(after):
        raise PoseReconstructionError(
            "Atom counts differ during coordinate comparison."
        )

    squared_distances: list[float] = []
    distances: list[float] = []

    for original, reconstructed in zip(before, after):
        dx = original[0] - reconstructed[0]
        dy = original[1] - reconstructed[1]
        dz = original[2] - reconstructed[2]

        squared_distance = dx * dx + dy * dy + dz * dz
        distance = math.sqrt(squared_distance)

        squared_distances.append(squared_distance)
        distances.append(distance)

    rmsd = math.sqrt(
        sum(squared_distances) / len(squared_distances)
    )

    return {
        "coordinate_rmsd_angstrom": round(rmsd, 8),
        "maximum_atom_displacement_angstrom": round(
            max(distances),
            8,
        ),
    }


def reconstruct_crystallographic_pose(
    template_path: Path,
    crystallographic_pdb_path: Path,
    output_sdf_path: Path,
) -> dict[str, Any]:
    """Transfer template bond orders onto crystallographic coordinates."""

    template = load_template_molecule(template_path)
    crystallographic_pose = load_crystallographic_pose(
        crystallographic_pdb_path
    )

    template_heavy_atoms = template.GetNumHeavyAtoms()
    pose_heavy_atoms = crystallographic_pose.GetNumHeavyAtoms()

    if template_heavy_atoms != pose_heavy_atoms:
        raise PoseReconstructionError(
            "Template and crystallographic pose have different "
            f"heavy-atom counts: template={template_heavy_atoms}, "
            f"pose={pose_heavy_atoms}."
        )

    original_coordinates = get_coordinates(
        crystallographic_pose
    )

    try:
        corrected_pose = AllChem.AssignBondOrdersFromTemplate(
            template,
            crystallographic_pose,
        )
    except Exception as error:
        raise PoseReconstructionError(
            "Could not assign template bond orders to the "
            "crystallographic ligand. Confirm that the template "
            "and extracted ligand represent the same compound."
        ) from error

    try:
        Chem.SanitizeMol(corrected_pose)
    except Exception as error:
        raise PoseReconstructionError(
            "The reconstructed ligand failed RDKit sanitization."
        ) from error

    corrected_coordinates = get_coordinates(corrected_pose)

    coordinate_check = calculate_coordinate_difference(
        before=original_coordinates,
        after=corrected_coordinates,
    )

    corrected_pose.SetProp("_Name", "1HSG_MK1_crystallographic")
    corrected_pose.SetProp(
        "DockForgeAI_Source",
        str(crystallographic_pdb_path),
    )
    corrected_pose.SetProp(
        "DockForgeAI_Bond_Order_Template",
        str(template_path),
    )

    output_sdf_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    writer = Chem.SDWriter(str(output_sdf_path))

    try:
        writer.write(corrected_pose)
    finally:
        writer.close()

    if (
        not output_sdf_path.exists()
        or output_sdf_path.stat().st_size == 0
    ):
        raise PoseReconstructionError(
            "The reconstructed SDF file was not created."
        )

    report = {
        "template_file": str(template_path),
        "crystallographic_source_file": str(
            crystallographic_pdb_path
        ),
        "output_file": str(output_sdf_path),
        "template_heavy_atoms": template_heavy_atoms,
        "crystallographic_heavy_atoms": pose_heavy_atoms,
        "output_heavy_atoms": corrected_pose.GetNumHeavyAtoms(),
        "output_total_atoms": corrected_pose.GetNumAtoms(),
        "output_bonds": corrected_pose.GetNumBonds(),
        "molecular_formula": (
            rdMolDescriptors.CalcMolFormula(corrected_pose)
        ),
        "formal_charge": Chem.GetFormalCharge(corrected_pose),
        "canonical_smiles": Chem.MolToSmiles(
            Chem.RemoveHs(corrected_pose),
            canonical=True,
        ),
        "sanitization_status": "passed",
        **coordinate_check,
        "coordinates_preserved": (
            coordinate_check[
                "maximum_atom_displacement_angstrom"
            ]
            < 1e-6
        ),
    }

    return report


def save_reconstruction_report(
    report: dict[str, Any],
    report_name: str,
) -> Path:
    """Save a pose-reconstruction report."""

    settings = load_settings()

    report_directory = (
        PROJECT_ROOT
        / settings["directories"]["ligand_reports"]
    )

    report_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    report_path = (
        report_directory
        / f"{report_name}_reconstruction.json"
    )

    with report_path.open("w", encoding="utf-8") as handle:
        json.dump(
            report,
            handle,
            indent=2,
        )

    return report_path


def reconstruct_mk1_pose() -> tuple[dict[str, Any], Path]:
    """Reconstruct the MK1 crystallographic pose from project files."""

    settings = load_settings()

    template_path = (
        PROJECT_ROOT
        / settings["directories"]["raw_ligands"]
        / "MK1_ideal.sdf"
    )

    crystallographic_path = (
        PROJECT_ROOT
        / settings["directories"]["extracted_ligands"]
        / "1HSG_MK1_chain-B_res-902.pdb"
    )

    output_path = (
        PROJECT_ROOT
        / settings["directories"]["crystallographic_ligands"]
        / "1HSG_MK1_crystallographic_corrected.sdf"
    )

    report = reconstruct_crystallographic_pose(
        template_path=template_path,
        crystallographic_pdb_path=crystallographic_path,
        output_sdf_path=output_path,
    )

    report_path = save_reconstruction_report(
        report=report,
        report_name="1HSG_MK1",
    )

    return report, report_path
