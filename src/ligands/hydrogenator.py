"""Add explicit hydrogens while preserving ligand heavy-atom coordinates."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import yaml
from rdkit import Chem
from rdkit.Chem import AllChem


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SETTINGS_PATH = PROJECT_ROOT / "config" / "settings.yaml"


class LigandHydrogenationError(RuntimeError):
    """Raised when ligand hydrogenation fails."""


def load_settings(
    settings_path: Path = SETTINGS_PATH,
) -> dict[str, Any]:
    """Load project configuration."""

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


def load_molecule(sdf_path: Path) -> Chem.Mol:
    """Load one sanitized molecule from an SDF file."""

    if not sdf_path.exists():
        raise FileNotFoundError(
            f"SDF file not found: {sdf_path}"
        )

    supplier = Chem.SDMolSupplier(
        str(sdf_path),
        removeHs=False,
        sanitize=True,
    )

    molecule = next(
        (
            mol
            for mol in supplier
            if mol is not None
        ),
        None,
    )

    if molecule is None:
        raise LigandHydrogenationError(
            f"Could not read molecule: {sdf_path}"
        )

    if molecule.GetNumConformers() == 0:
        raise LigandHydrogenationError(
            "The molecule has no 3D conformer."
        )

    return molecule


def get_heavy_atom_coordinates(
    molecule: Chem.Mol,
) -> list[tuple[float, float, float]]:
    """Return heavy-atom coordinates in atom order."""

    conformer = molecule.GetConformer()
    coordinates: list[tuple[float, float, float]] = []

    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue

        position = conformer.GetAtomPosition(atom.GetIdx())

        coordinates.append(
            (
                float(position.x),
                float(position.y),
                float(position.z),
            )
        )

    return coordinates


def calculate_coordinate_change(
    before: list[tuple[float, float, float]],
    after: list[tuple[float, float, float]],
) -> dict[str, float]:
    """Calculate heavy-atom RMSD and maximum displacement."""

    if len(before) != len(after):
        raise LigandHydrogenationError(
            "Heavy-atom counts differ before and after hydrogenation."
        )

    squared_distances: list[float] = []
    distances: list[float] = []

    for original, final in zip(before, after):
        dx = original[0] - final[0]
        dy = original[1] - final[1]
        dz = original[2] - final[2]

        squared_distance = dx * dx + dy * dy + dz * dz

        squared_distances.append(squared_distance)
        distances.append(math.sqrt(squared_distance))

    rmsd = math.sqrt(
        sum(squared_distances) / len(squared_distances)
    )

    return {
        "heavy_atom_rmsd_angstrom": round(rmsd, 8),
        "maximum_heavy_atom_displacement_angstrom": round(
            max(distances),
            8,
        ),
    }


def add_hydrogens_preserving_pose(
    input_sdf: Path,
    output_sdf: Path,
) -> dict[str, Any]:
    """Add explicit hydrogens and optimize only hydrogen coordinates."""

    molecule = load_molecule(input_sdf)
    molecule = Chem.RemoveHs(molecule)

    original_coordinates = get_heavy_atom_coordinates(molecule)

    hydrogenated = Chem.AddHs(
        molecule,
        addCoords=True,
    )

    properties = AllChem.MMFFGetMoleculeProperties(
        hydrogenated,
        mmffVariant="MMFF94",
    )

    if properties is None:
        raise LigandHydrogenationError(
            "MMFF94 parameters are unavailable for this ligand."
        )

    force_field = AllChem.MMFFGetMoleculeForceField(
        hydrogenated,
        properties,
    )

    if force_field is None:
        raise LigandHydrogenationError(
            "Could not create the MMFF94 force field."
        )

    # Keep every heavy atom fixed.
    for atom in hydrogenated.GetAtoms():
        if atom.GetAtomicNum() != 1:
            force_field.AddFixedPoint(atom.GetIdx())

    result = force_field.Minimize(maxIts=500)

    optimization_status = (
        "converged"
        if result == 0
        else "maximum iterations reached"
    )

    final_coordinates = get_heavy_atom_coordinates(hydrogenated)

    coordinate_change = calculate_coordinate_change(
        before=original_coordinates,
        after=final_coordinates,
    )

    Chem.SanitizeMol(hydrogenated)

    hydrogenated.SetProp(
        "_Name",
        "1HSG_MK1_crystallographic_H",
    )

    output_sdf.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    writer = Chem.SDWriter(str(output_sdf))

    try:
        writer.write(hydrogenated)
    finally:
        writer.close()

    if not output_sdf.exists() or output_sdf.stat().st_size == 0:
        raise LigandHydrogenationError(
            "Hydrogenated ligand was not written."
        )

    hydrogen_count = sum(
        1
        for atom in hydrogenated.GetAtoms()
        if atom.GetAtomicNum() == 1
    )

    return {
        "input_file": str(input_sdf),
        "output_file": str(output_sdf),
        "heavy_atoms": hydrogenated.GetNumHeavyAtoms(),
        "total_atoms": hydrogenated.GetNumAtoms(),
        "explicit_hydrogens": hydrogen_count,
        "formal_charge": Chem.GetFormalCharge(hydrogenated),
        "hydrogen_optimization": optimization_status,
        **coordinate_change,
        "heavy_atom_coordinates_preserved": (
            coordinate_change[
                "maximum_heavy_atom_displacement_angstrom"
            ]
            < 1e-6
        ),
    }


def prepare_mk1_hydrogens() -> tuple[dict[str, Any], Path]:
    """Prepare hydrogenated crystallographic MK1."""

    settings = load_settings()

    input_path = (
        PROJECT_ROOT
        / settings["directories"]["crystallographic_ligands"]
        / "1HSG_MK1_crystallographic_corrected.sdf"
    )

    output_path = (
        PROJECT_ROOT
        / settings["directories"]["docking_ready_ligands"]
        / "1HSG_MK1_crystallographic_H.sdf"
    )

    report = add_hydrogens_preserving_pose(
        input_sdf=input_path,
        output_sdf=output_path,
    )

    report_directory = (
        PROJECT_ROOT
        / settings["directories"]["docking_reports"]
    )

    report_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    report_path = (
        report_directory
        / "1HSG_MK1_hydrogenation.json"
    )

    with report_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)

    return report, report_path
