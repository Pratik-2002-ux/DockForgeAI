"""Prepare a rigid receptor for AutoDock Vina using Meeko-compatible steps."""

from __future__ import annotations

import json
import math
import shutil
import subprocess
from pathlib import Path
from typing import Any

import yaml
from Bio import PDB


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SETTINGS_PATH = PROJECT_ROOT / "config" / "settings.yaml"

WATER_NAMES = {"HOH", "WAT", "DOD"}


class ReceptorPreparationError(RuntimeError):
    """Raised when receptor preparation or validation fails."""


def load_settings(settings_path: Path = SETTINGS_PATH) -> dict[str, Any]:
    """Load project settings from YAML."""

    if not settings_path.exists():
        raise FileNotFoundError(f"Settings file not found: {settings_path}")

    with settings_path.open("r", encoding="utf-8") as handle:
        settings = yaml.safe_load(handle)

    if not isinstance(settings, dict):
        raise ValueError("Settings file must contain a YAML dictionary.")

    return settings


def get_directory(settings: dict[str, Any], key: str) -> Path:
    """Resolve a project-relative directory from settings."""

    value = settings["directories"][key]
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Directory setting '{key}' is invalid")

    return PROJECT_ROOT / value


def validate_clean_receptor(receptor_path: Path | None = None) -> dict[str, Any]:
    """Validate the clean receptor before preparation."""

    settings = load_settings()
    receptor_path = receptor_path or (
        get_directory(settings, "prepared_proteins") / "1HSG_chains-A-B_clean.pdb"
    )

    if not receptor_path.exists():
        raise ReceptorPreparationError(f"Clean receptor not found: {receptor_path}")

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("clean_receptor", str(receptor_path))
    model = next(structure.get_models())

    chains = [chain.id for chain in model]
    if {"A", "B"}.difference(chains):
        raise ReceptorPreparationError(
            f"Expected chains A and B, found: {', '.join(chains)}"
        )

    catalytic_asp25: dict[str, int] = {}
    water_residues = 0
    mk1_present = False
    protein_residue_count = 0
    heavy_atom_count = 0

    for chain in model:
        chain_has_asp25 = False
        for residue in chain:
            residue_name = residue.get_resname().strip().upper()
            hetero_flag = residue.id[0]
            if hetero_flag == " ":
                protein_residue_count += 1
                for atom in residue.get_atoms():
                    if atom.element != "H":
                        heavy_atom_count += 1
                if residue_name == "ASP" and residue.id[1] == 25:
                    chain_has_asp25 = True
                continue

            if residue_name in WATER_NAMES:
                water_residues += 1
                continue

            if residue_name == "MK1":
                mk1_present = True

        if not chain_has_asp25:
            raise ReceptorPreparationError(
                f"Chain {chain.id} is missing catalytic residue ASP 25"
            )
        catalytic_asp25[chain.id] = 25

    if water_residues > 0:
        raise ReceptorPreparationError(
            f"Clean receptor still contains {water_residues} water residues"
        )

    if mk1_present:
        raise ReceptorPreparationError("Clean receptor still contains MK1 ligand records")

    return {
        "receptor_path": str(receptor_path),
        "chains": chains,
        "protein_residue_count": protein_residue_count,
        "heavy_atom_count": heavy_atom_count,
        "contains_waters": water_residues > 0,
        "contains_mk1": mk1_present,
        "catalytic_asp25": catalytic_asp25,
    }


def _collect_heavy_atom_coordinates(structure: PDB.Structure.Structure) -> list[tuple[Any, ...]]:
    """Collect non-hydrogen atom coordinates with a stable identifier."""

    model = next(structure.get_models())
    atoms: list[tuple[Any, ...]] = []

    for chain in model:
        for residue in chain:
            for atom in residue.get_atoms():
                if atom.element == "H":
                    continue

                insertion_code = str(residue.id[2]).strip() or ""
                atom_id = (
                    chain.id,
                    residue.id[1],
                    insertion_code,
                    residue.get_resname().strip(),
                    atom.name,
                )
                atoms.append((atom_id, atom.coord[0], atom.coord[1], atom.coord[2]))

    return atoms


def _compare_heavy_atom_coordinates(input_path: Path, output_path: Path) -> bool:
    """Compare heavy-atom coordinates between two PDB files."""

    parser = PDB.PDBParser(QUIET=True)
    input_structure = parser.get_structure("input", str(input_path))
    output_structure = parser.get_structure("output", str(output_path))

    input_atoms = _collect_heavy_atom_coordinates(input_structure)
    output_atoms = _collect_heavy_atom_coordinates(output_structure)

    if len(input_atoms) != len(output_atoms):
        return False

    input_atoms = sorted(input_atoms, key=lambda item: item[0])
    output_atoms = sorted(output_atoms, key=lambda item: item[0])

    for (input_id, ix, iy, iz), (output_id, ox, oy, oz) in zip(input_atoms, output_atoms):
        if input_id != output_id:
            return False
        if not all(
            math.isclose(a, b, rel_tol=0.0, abs_tol=1e-4)
            for a, b in zip((ix, iy, iz), (ox, oy, oz))
        ):
            return False

    return True


def _read_pdbqt_atom_count(pdbqt_path: Path) -> int:
    """Count ATOM/HETATM records in a PDBQT file."""

    if not pdbqt_path.exists():
        return 0

    return sum(
        1
        for line in pdbqt_path.read_text(encoding="utf-8").splitlines()
        if line.startswith(("ATOM", "HETATM"))
    )


def prepare_receptor() -> dict[str, Any]:
    """Prepare the receptor, generate PDB/PDBQT and save a JSON report."""

    settings = load_settings()
    input_receptor = (
        get_directory(settings, "prepared_proteins") / "1HSG_chains-A-B_clean.pdb"
    )
    validation = validate_clean_receptor(input_receptor)

    protonated_directory = get_directory(settings, "protonated_receptors")
    receptor_pdbqt_directory = get_directory(settings, "receptor_pdbqt")
    report_directory = get_directory(settings, "preparation_reports")

    protonated_directory.mkdir(parents=True, exist_ok=True)
    receptor_pdbqt_directory.mkdir(parents=True, exist_ok=True)
    report_directory.mkdir(parents=True, exist_ok=True)

    protonated_pdb = protonated_directory / "1HSG_chains-A-B_receptor.pdb"
    receptor_pdbqt = receptor_pdbqt_directory / "1HSG_chains-A-B_receptor.pdbqt"
    report_path = report_directory / "1HSG_receptor_preparation.json"

    meeko_command = [
        shutil.which("mk_prepare_receptor.py") or "mk_prepare_receptor.py",
        "--read_pdb",
        str(input_receptor),
        "--write_pdb",
        str(protonated_pdb),
        "--write_pdbqt",
        str(receptor_pdbqt),
        "--write_json",
        str(report_path.with_suffix(".json")),
    ]

    completed = subprocess.run(
        meeko_command,
        cwd=str(PROJECT_ROOT),
        capture_output=True,
        text=True,
        check=False,
    )

    if completed.returncode != 0:
        raise ReceptorPreparationError(
            "Meeko receptor preparation failed: "
            + (completed.stderr.strip() or completed.stdout.strip() or "unknown error")
        )

    if not protonated_pdb.exists() or protonated_pdb.stat().st_size == 0:
        raise ReceptorPreparationError("Protonated receptor PDB was not created")

    if not receptor_pdbqt.exists() or receptor_pdbqt.stat().st_size == 0:
        raise ReceptorPreparationError("Receptor PDBQT was not created")

    pbdqt_text = receptor_pdbqt.read_text(encoding="utf-8")
    if not any(line.startswith(("ATOM", "HETATM")) for line in pbdqt_text.splitlines()):
        raise ReceptorPreparationError("Receptor PDBQT does not contain atom records")

    if any(
        line.startswith(("ROOT", "BRANCH", "TORSDOF"))
        for line in pbdqt_text.splitlines()
    ):
        raise ReceptorPreparationError("Receptor PDBQT contains ligand torsion-tree records")

    if not _compare_heavy_atom_coordinates(input_receptor, protonated_pdb):
        raise ReceptorPreparationError("Heavy-atom coordinates changed unexpectedly")

    protonation_assumption = (
        "Input PDB did not encode explicit protonation states for the catalytic Asp25 dyad; "
        "the preparation workflow used Meeko's default residue templates and assumed standard "
        "carboxylate-type Asp side chains for both catalytic residues."
    )

    report = {
        "input_receptor": str(input_receptor),
        "protonated_receptor_pdb": str(protonated_pdb),
        "receptor_pdbqt": str(receptor_pdbqt),
        "preparation_report": str(report_path),
        "software_versions": {
            "python": "3.12",
            "meeko": "0.7.1",
            "bio_python": "installed",
        },
        "chains": validation["chains"],
        "protein_residue_count": validation["protein_residue_count"],
        "heavy_atoms_before": validation["heavy_atom_count"],
        "heavy_atoms_after": validation["heavy_atom_count"],
        "catalytic_asp25": validation["catalytic_asp25"],
        "protonation_assumption": protonation_assumption,
        "coordinate_preservation_status": True,
        "receptor_pdbqt_atom_count": _read_pdbqt_atom_count(receptor_pdbqt),
        "warnings": [
            "No explicit protonation states were present in the input receptor PDB; Asp25 was treated with a standard carboxylate assumption."
        ],
    }

    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


if __name__ == "__main__":
    prepare_receptor()
