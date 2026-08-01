"""Export AutoDock Vina PDBQT poses to a chemically valid SDF file.

The exporter reconstructs each ligand from Meeko's embedded SMILES and
SMILES-to-PDBQT atom-index mapping. This preserves the chemical graph while
transferring the actual coordinates of every Vina pose.

It intentionally does not infer bond orders from PDBQT coordinates.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import yaml
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Geometry import Point3D


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SETTINGS_PATH = PROJECT_ROOT / "config" / "settings.yaml"


class PoseExportError(RuntimeError):
    """Raised when docked poses cannot be exported safely."""


def load_settings(
    settings_path: Path = SETTINGS_PATH,
) -> dict[str, Any]:
    """Load project settings from YAML."""

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


def _load_reference_molecule(reference_sdf: Path) -> Chem.Mol:
    """Load the chemically corrected crystallographic reference ligand."""

    if not reference_sdf.exists():
        raise FileNotFoundError(
            f"Reference SDF not found: {reference_sdf}"
        )

    supplier = Chem.SDMolSupplier(
        str(reference_sdf),
        removeHs=True,
        sanitize=True,
    )

    molecules = [
        molecule
        for molecule in supplier
        if molecule is not None
    ]

    if len(molecules) != 1:
        raise PoseExportError(
            "Reference SDF must contain exactly one valid molecule."
        )

    reference = Chem.Mol(molecules[0])

    if reference.GetNumConformers() == 0:
        raise PoseExportError(
            "Reference SDF does not contain coordinates."
        )

    return reference


def _parse_meeko_smiles(lines: list[str]) -> str:
    """Extract the ligand SMILES embedded by Meeko in a PDBQT file."""

    smiles_parts: list[str] = []

    for line in lines:
        if line.startswith("REMARK SMILES IDX"):
            continue

        if line.startswith("REMARK SMILES "):
            value = line[len("REMARK SMILES "):].strip()

            if value:
                smiles_parts.append(value)

    if not smiles_parts:
        raise PoseExportError(
            "PDBQT does not contain a Meeko REMARK SMILES record."
        )

    return "".join(smiles_parts)


def _parse_smiles_index_mapping(
    lines: list[str],
    molecule_atom_count: int,
) -> dict[int, int]:
    """Parse Meeko SMILES atom indices and PDBQT atom serial numbers.

    Returns
    -------
    dict
        Mapping from zero-based RDKit atom index to PDBQT atom serial number.
    """

    values: list[int] = []

    prefix = "REMARK SMILES IDX"

    for line in lines:
        if not line.startswith(prefix):
            continue

        tokens = line[len(prefix):].split()

        try:
            values.extend(int(token) for token in tokens)
        except ValueError as error:
            raise PoseExportError(
                f"Invalid REMARK SMILES IDX record: {line}"
            ) from error

    if not values:
        raise PoseExportError(
            "PDBQT does not contain REMARK SMILES IDX mapping."
        )

    if len(values) % 2 != 0:
        raise PoseExportError(
            "REMARK SMILES IDX contains an odd number of indices."
        )

    pairs = list(zip(values[0::2], values[1::2]))

    expected_smiles_indices = set(
        range(1, molecule_atom_count + 1)
    )

    first_values = {first for first, _ in pairs}
    second_values = {second for _, second in pairs}

    if expected_smiles_indices.issubset(first_values):
        mapping = {
            smiles_index - 1: pdbqt_serial
            for smiles_index, pdbqt_serial in pairs
            if smiles_index <= molecule_atom_count
        }

    elif expected_smiles_indices.issubset(second_values):
        # Defensive support in case a Meeko release reverses pair order.
        mapping = {
            smiles_index - 1: pdbqt_serial
            for pdbqt_serial, smiles_index in pairs
            if smiles_index <= molecule_atom_count
        }

    else:
        raise PoseExportError(
            "Could not determine the direction of the SMILES atom mapping. "
            f"Expected atom indices 1–{molecule_atom_count}."
        )

    if set(mapping) != set(range(molecule_atom_count)):
        missing = sorted(
            set(range(molecule_atom_count)) - set(mapping)
        )

        raise PoseExportError(
            "SMILES atom mapping is incomplete. Missing zero-based "
            f"atom indices: {missing}"
        )

    if len(set(mapping.values())) != len(mapping):
        raise PoseExportError(
            "Multiple SMILES atoms map to the same PDBQT atom."
        )

    return mapping


def _parse_pdbqt_coordinates(
    lines: list[str],
) -> dict[int, tuple[float, float, float]]:
    """Extract PDBQT coordinates keyed by atom serial number."""

    coordinates: dict[int, tuple[float, float, float]] = {}

    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue

        try:
            serial = int(line[6:11].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
        except ValueError:
            # Defensive fallback for non-standard spacing.
            parts = line.split()

            if len(parts) < 8:
                raise PoseExportError(
                    f"Could not parse PDBQT atom line: {line}"
                )

            try:
                serial = int(parts[1])
                x = float(parts[5])
                y = float(parts[6])
                z = float(parts[7])
            except ValueError as error:
                raise PoseExportError(
                    f"Could not parse PDBQT atom line: {line}"
                ) from error

        coordinates[serial] = (x, y, z)

    if not coordinates:
        raise PoseExportError(
            "PDBQT pose contains no coordinate records."
        )

    return coordinates


def _build_pose_molecule(
    pose_path: Path,
    reference_molecule: Chem.Mol,
    pose_number: int,
) -> Chem.Mol:
    """Reconstruct one chemically valid pose from a PDBQT file."""

    if not pose_path.exists():
        raise FileNotFoundError(
            f"Pose PDBQT not found: {pose_path}"
        )

    lines = pose_path.read_text(
        encoding="utf-8",
    ).splitlines()

    smiles = _parse_meeko_smiles(lines)

    molecule = Chem.MolFromSmiles(
        smiles,
        sanitize=True,
    )

    if molecule is None:
        raise PoseExportError(
            f"RDKit could not parse Meeko SMILES from {pose_path}."
        )

    molecule = Chem.RemoveHs(molecule)

    if molecule.GetNumAtoms() != reference_molecule.GetNumHeavyAtoms():
        raise PoseExportError(
            f"Pose {pose_number} SMILES has {molecule.GetNumAtoms()} "
            "heavy atoms; expected "
            f"{reference_molecule.GetNumHeavyAtoms()}."
        )

    mapping = _parse_smiles_index_mapping(
        lines=lines,
        molecule_atom_count=molecule.GetNumAtoms(),
    )

    pdbqt_coordinates = _parse_pdbqt_coordinates(lines)

    # Remove every conformer before assigning the docked coordinates.
    molecule.RemoveAllConformers()

    conformer = Chem.Conformer(molecule.GetNumAtoms())
    conformer.Set3D(True)

    for rdkit_atom_index in range(molecule.GetNumAtoms()):
        pdbqt_serial = mapping[rdkit_atom_index]

        if pdbqt_serial not in pdbqt_coordinates:
            raise PoseExportError(
                f"Pose {pose_number}: PDBQT atom serial "
                f"{pdbqt_serial} is missing."
            )

        x, y, z = pdbqt_coordinates[pdbqt_serial]

        conformer.SetAtomPosition(
            rdkit_atom_index,
            Point3D(x, y, z),
        )

    molecule.AddConformer(
        conformer,
        assignId=True,
    )

    Chem.AssignStereochemistry(
        molecule,
        cleanIt=True,
        force=True,
    )
    Chem.SanitizeMol(molecule)

    molecule.SetProp(
        "_Name",
        f"1HSG_MK1_pose_{pose_number:03d}",
    )
    molecule.SetIntProp(
        "pose_rank",
        pose_number,
    )
    molecule.SetProp(
        "source_pdbqt",
        str(pose_path),
    )
    molecule.SetProp(
        "export_method",
        "rdkit_meeko_smiles_idx",
    )

    reference_formula = rdMolDescriptors.CalcMolFormula(
        reference_molecule
    )
    pose_formula = rdMolDescriptors.CalcMolFormula(
        molecule
    )

    if pose_formula != reference_formula:
        raise PoseExportError(
            f"Pose {pose_number} formula is {pose_formula}; "
            f"expected {reference_formula}."
        )

    reference_smiles = Chem.MolToSmiles(
        reference_molecule,
        canonical=True,
        isomericSmiles=True,
    )
    pose_smiles = Chem.MolToSmiles(
        molecule,
        canonical=True,
        isomericSmiles=True,
    )

    if pose_smiles != reference_smiles:
        # Some source SDFs may differ only in stereochemical annotation.
        reference_connectivity = Chem.MolToSmiles(
            reference_molecule,
            canonical=True,
            isomericSmiles=False,
        )
        pose_connectivity = Chem.MolToSmiles(
            molecule,
            canonical=True,
            isomericSmiles=False,
        )

        if pose_connectivity != reference_connectivity:
            raise PoseExportError(
                f"Pose {pose_number} chemical graph does not match "
                "the crystallographic reference."
            )

    return molecule


def _coordinate_signature(
    molecule: Chem.Mol,
    precision: int = 4,
) -> tuple[tuple[float, float, float], ...]:
    """Return a rounded coordinate signature for duplicate detection."""

    if molecule.GetNumConformers() != 1:
        raise PoseExportError(
            "Pose molecule must have exactly one conformer."
        )

    conformer = molecule.GetConformer()

    return tuple(
        (
            round(conformer.GetAtomPosition(index).x, precision),
            round(conformer.GetAtomPosition(index).y, precision),
            round(conformer.GetAtomPosition(index).z, precision),
        )
        for index in range(molecule.GetNumAtoms())
    )


def _coordinate_rmsd_by_index(
    first: Chem.Mol,
    second: Chem.Mol,
) -> float:
    """Calculate direct coordinate RMSD using the same atom order."""

    if first.GetNumAtoms() != second.GetNumAtoms():
        raise PoseExportError(
            "Cannot compare poses with different atom counts."
        )

    first_conf = first.GetConformer()
    second_conf = second.GetConformer()

    squared_sum = 0.0

    for atom_index in range(first.GetNumAtoms()):
        first_position = first_conf.GetAtomPosition(atom_index)
        second_position = second_conf.GetAtomPosition(atom_index)

        dx = first_position.x - second_position.x
        dy = first_position.y - second_position.y
        dz = first_position.z - second_position.z

        squared_sum += dx * dx + dy * dy + dz * dz

    return math.sqrt(
        squared_sum / first.GetNumAtoms()
    )


def export_poses(
    pose_pdbqts: list[Path],
    output_sdf: Path,
    reference_sdf: Path,
    report_path: Path,
) -> dict[str, Any]:
    """Export all docked poses to one multi-molecule SDF."""

    if not pose_pdbqts:
        raise PoseExportError(
            "No PDBQT poses were supplied for export."
        )

    reference_molecule = _load_reference_molecule(
        reference_sdf
    )

    output_sdf.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    report_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    exported_molecules: list[Chem.Mol] = []

    for pose_number, pose_path in enumerate(
        pose_pdbqts,
        start=1,
    ):
        exported_molecules.append(
            _build_pose_molecule(
                pose_path=pose_path,
                reference_molecule=reference_molecule,
                pose_number=pose_number,
            )
        )

    signatures = {
        _coordinate_signature(molecule)
        for molecule in exported_molecules
    }

    if len(exported_molecules) > 1 and len(signatures) == 1:
        raise PoseExportError(
            "All exported poses have identical coordinates. "
            "Pose coordinate transfer failed."
        )

    pose_1_rmsds = [
        round(
            _coordinate_rmsd_by_index(
                exported_molecules[0],
                molecule,
            ),
            6,
        )
        for molecule in exported_molecules
    ]

    temporary_path = output_sdf.with_suffix(
        output_sdf.suffix + ".tmp"
    )

    writer = Chem.SDWriter(str(temporary_path))

    try:
        for molecule in exported_molecules:
            writer.write(molecule)
    finally:
        writer.close()

    if (
        not temporary_path.exists()
        or temporary_path.stat().st_size == 0
    ):
        raise PoseExportError(
            "Temporary multi-pose SDF was not created."
        )

    # Validate the written file before replacing the final output.
    supplier = Chem.SDMolSupplier(
        str(temporary_path),
        removeHs=True,
        sanitize=True,
    )

    written_molecules = [
        molecule
        for molecule in supplier
        if molecule is not None
    ]

    if len(written_molecules) != len(exported_molecules):
        temporary_path.unlink(missing_ok=True)

        raise PoseExportError(
            "Written SDF pose count does not match input pose count."
        )

    temporary_path.replace(output_sdf)

    report = {
        "export_method": "rdkit_meeko_smiles_idx",
        "output_sdf": str(output_sdf),
        "pose_count": len(exported_molecules),
        "reference_sdf": str(reference_sdf),
        "unique_coordinate_sets": len(signatures),
        "pose_rmsd_relative_to_pose_1_angstrom": pose_1_rmsds,
        "all_poses_identical": len(signatures) == 1,
    }

    report_path.write_text(
        json.dumps(report, indent=2),
        encoding="utf-8",
    )

    return report
