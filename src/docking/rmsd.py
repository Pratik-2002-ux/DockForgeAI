"""Calculate crystallographic RMSD for each docked pose."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import yaml
from rdkit import Chem
from rdkit.Chem import rdMolAlign, rdMolDescriptors


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SETTINGS_PATH = PROJECT_ROOT / "config" / "settings.yaml"


class RMSDCalculationError(RuntimeError):
    """Raised when RMSD calculation fails."""


def load_settings(settings_path: Path = SETTINGS_PATH) -> dict[str, Any]:
    """Load project settings from YAML."""

    if not settings_path.exists():
        raise FileNotFoundError(f"Settings file not found: {settings_path}")

    with settings_path.open("r", encoding="utf-8") as handle:
        settings = yaml.safe_load(handle)

    if not isinstance(settings, dict):
        raise ValueError("Settings file must contain a YAML dictionary.")

    return settings


def _prepare_molecule(molecule: Chem.Mol) -> Chem.Mol:
    """Return a copy of a molecule without explicit hydrogens and with a single conformer."""

    prepared = Chem.Mol(molecule)
    prepared = Chem.RemoveHs(prepared)
    if prepared.GetNumConformers() == 0:
        raise RMSDCalculationError("Molecule has no conformer")
    return prepared


def calculate_rmsd(
    reference_sdf: Path,
    pose_sdf: Path,
    score_table_path: Path,
    output_table_path: Path,
    metadata_path: Path,
    threshold: float,
) -> dict[str, Any]:
    """Merge Vina scores with crystallographic RMSD values."""

    settings = load_settings()
    reference_supplier = Chem.SDMolSupplier(str(reference_sdf), removeHs=False, sanitize=True)
    reference_molecules = [mol for mol in reference_supplier if mol is not None]
    if len(reference_molecules) != 1:
        raise RMSDCalculationError("Reference SDF should contain exactly one molecule")
    reference_molecule = _prepare_molecule(reference_molecules[0])

    pose_supplier = Chem.SDMolSupplier(str(pose_sdf), removeHs=False, sanitize=True)
    pose_molecules = [mol for mol in pose_supplier if mol is not None]
    if len(pose_molecules) == 0:
        raise RMSDCalculationError("Exported pose SDF contains no molecules")

    with score_table_path.open("r", encoding="utf-8", newline="") as handle:
        score_rows = list(csv.DictReader(handle))

    if len(score_rows) != len(pose_molecules):
        raise RMSDCalculationError("Pose count mismatch between score rows and exported molecules")

    results: list[dict[str, Any]] = []
    for row, probe in zip(score_rows, pose_molecules):
        probe_molecule = _prepare_molecule(probe)
        if probe_molecule.GetNumHeavyAtoms() != reference_molecule.GetNumHeavyAtoms():
            raise RMSDCalculationError("Pose and reference heavy-atom counts differ")
        if rdMolDescriptors.CalcMolFormula(probe_molecule) != rdMolDescriptors.CalcMolFormula(reference_molecule):
            raise RMSDCalculationError("Pose and reference formulas differ")
        if Chem.GetFormalCharge(probe_molecule) != Chem.GetFormalCharge(reference_molecule):
            raise RMSDCalculationError("Pose and reference formal charges differ")

        try:
            in_place_rmsd = rdMolAlign.CalcRMS(reference_molecule, probe_molecule)
        except Exception as error:
            raise RMSDCalculationError(f"RMSD calculation failed: {error}") from error

        try:
            aligned_best_rmsd = rdMolAlign.GetBestRMS(reference_molecule, probe_molecule)
        except Exception as error:
            aligned_best_rmsd = float("nan")

        in_place_rmsd = float(in_place_rmsd)
        aligned_best_rmsd = float(aligned_best_rmsd)
        if not (in_place_rmsd >= 0.0 and in_place_rmsd == in_place_rmsd):
            raise RMSDCalculationError("Calculated RMSD is not finite and non-negative")

        row.update(
            {
                "crystal_in_place_rmsd_angstrom": round(in_place_rmsd, 3),
                "aligned_best_rmsd_angstrom": round(aligned_best_rmsd, 3),
                "passes_2A_threshold": in_place_rmsd <= threshold,
            }
        )
        results.append(row)

    output_table_path.parent.mkdir(parents=True, exist_ok=True)
    with output_table_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "pose_rank",
                "affinity_kcal_mol",
                "vina_rmsd_lower_bound_angstrom",
                "vina_rmsd_upper_bound_angstrom",
                "crystal_in_place_rmsd_angstrom",
                "aligned_best_rmsd_angstrom",
                "passes_2A_threshold",
                "atom_records",
                "pose_pdbqt",
            ],
        )
        writer.writeheader()
        writer.writerows(results)

    metadata = {
        "reference_sdf": str(reference_sdf),
        "pose_sdf": str(pose_sdf),
        "pose_count": len(results),
        "threshold_angstrom": threshold,
    }
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return metadata
