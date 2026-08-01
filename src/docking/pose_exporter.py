"""Export docked Vina poses to an SDF file."""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path
from typing import Any

import yaml
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SETTINGS_PATH = PROJECT_ROOT / "config" / "settings.yaml"


class PoseExportError(RuntimeError):
    """Raised when docked poses cannot be exported to SDF."""


def load_settings(settings_path: Path = SETTINGS_PATH) -> dict[str, Any]:
    """Load project settings from YAML."""

    if not settings_path.exists():
        raise FileNotFoundError(f"Settings file not found: {settings_path}")

    with settings_path.open("r", encoding="utf-8") as handle:
        settings = yaml.safe_load(handle)

    if not isinstance(settings, dict):
        raise ValueError("Settings file must contain a YAML dictionary.")

    return settings


def _load_reference_molecule(reference_sdf: Path) -> Chem.Mol:
    """Load the chemically corrected reference ligand SDF."""

    supplier = Chem.SDMolSupplier(str(reference_sdf), removeHs=False, sanitize=True)
    molecules = [mol for mol in supplier if mol is not None]
    if len(molecules) != 1:
        raise PoseExportError("Reference SDF should contain exactly one molecule")
    reference = molecules[0]
    if reference.GetNumConformers() == 0:
        raise PoseExportError("Reference SDF has no conformer")
    return reference


def _fallback_export_from_pdbqt(
    pose_pdbqts: list[Path],
    reference_molecule: Chem.Mol,
    output_sdf: Path,
) -> list[Chem.Mol]:
    """Create an SDF from PDBQT poses using the reference template as a fallback."""

    exported: list[Chem.Mol] = []
    for pose_path in pose_pdbqts:
        lines = pose_path.read_text(encoding="utf-8").splitlines()
        atom_lines = [line for line in lines if line.startswith(("ATOM", "HETATM"))]
        if not atom_lines:
            raise PoseExportError(f"PDBQT pose lacks atom records: {pose_path}")

        template = Chem.Mol(reference_molecule)
        template = Chem.RemoveHs(template)
        heavy_atom_lines = [line for line in atom_lines if line.split()[0] in {"ATOM", "HETATM"} and line.split()[2] != "H"]
        if template.GetNumHeavyAtoms() != len(heavy_atom_lines):
            raise PoseExportError("Heavy-atom count mismatch between template and pose")

        conformer = Chem.Conformer(template.GetNumAtoms())
        for atom_index, line in enumerate(heavy_atom_lines):
            parts = line.split()
            if len(parts) < 8:
                raise PoseExportError(f"Unexpected atom line in {pose_path}: {line}")
            x = float(parts[5])
            y = float(parts[6])
            z = float(parts[7])
            conformer.SetAtomPosition(atom_index, Chem.rdGeometry.Point3D(x, y, z))
        template.AddConformer(conformer, assignId=True)
        template.SetProp("_Name", pose_path.stem)
        template.SetProp("export_method", "fallback_rdkit")
        exported.append(template)

    writer = Chem.SDWriter(str(output_sdf))
    try:
        for molecule in exported:
            writer.write(molecule)
    finally:
        writer.close()

    return exported


def export_poses(
    pose_pdbqts: list[Path],
    output_sdf: Path,
    reference_sdf: Path,
    report_path: Path,
) -> dict[str, Any]:
    """Export all poses to a single SDF and write a JSON report."""

    settings = load_settings()
    output_sdf.parent.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)

    meeko_executable = settings.get("docking", {}).get("meeko_export_executable", "mk_export.py")
    meeko_path = shutil.which(meeko_executable)
    if meeko_path is None:
        meeko_path = meeko_executable

    reference_molecule = _load_reference_molecule(reference_sdf)
    export_method = "meeko"

    if meeko_path:
        command = [meeko_path, "-s", str(output_sdf)]
        command.extend(str(path) for path in pose_pdbqts)
        completed = subprocess.run(
            command,
            cwd=str(PROJECT_ROOT),
            capture_output=True,
            text=True,
            check=False,
        )
        if completed.returncode == 0 and output_sdf.exists() and output_sdf.stat().st_size > 0:
            export_method = "meeko"
        else:
            export_method = "fallback_rdkit"
            _fallback_export_from_pdbqt(pose_pdbqts, reference_molecule, output_sdf)
    else:
        export_method = "fallback_rdkit"
        _fallback_export_from_pdbqt(pose_pdbqts, reference_molecule, output_sdf)

    supplier = Chem.SDMolSupplier(str(output_sdf), removeHs=False, sanitize=True)
    exported_molecules = [mol for mol in supplier if mol is not None]
    if len(exported_molecules) != len(pose_pdbqts):
        raise PoseExportError("Exported pose count does not match the number of input poses")

    for index, molecule in enumerate(exported_molecules, start=1):
        if molecule.GetNumConformers() != 1:
            raise PoseExportError(f"Pose {index} does not have a single conformer")
        if molecule.GetNumHeavyAtoms() != 45:
            raise PoseExportError(f"Pose {index} does not have 45 heavy atoms")
        if Chem.GetFormalCharge(molecule) != Chem.GetFormalCharge(reference_molecule):
            raise PoseExportError(f"Pose {index} has a different formal charge")
        formula = rdMolDescriptors.CalcMolFormula(molecule)
        reference_formula = rdMolDescriptors.CalcMolFormula(reference_molecule)
        if formula != reference_formula:
            raise PoseExportError(f"Pose {index} has formula {formula}, expected {reference_formula}")

    report = {
        "export_method": export_method,
        "output_sdf": str(output_sdf),
        "pose_count": len(exported_molecules),
        "reference_sdf": str(reference_sdf),
    }
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report
