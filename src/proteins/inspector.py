"""Protein structure inspection utilities for DockForgeAI.

This module examines a downloaded PDB structure before receptor preparation.
It reports chains, atoms, residues, waters, ligands, metal ions and selected
experimental metadata.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any

import yaml
from Bio.PDB import PDBParser


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SETTINGS_PATH = PROJECT_ROOT / "config" / "settings.yaml"


WATER_NAMES = {
    "HOH",
    "WAT",
    "DOD",
}

COMMON_METALS = {
    "LI",
    "NA",
    "K",
    "RB",
    "CS",
    "MG",
    "CA",
    "SR",
    "BA",
    "MN",
    "FE",
    "CO",
    "NI",
    "CU",
    "ZN",
    "CD",
    "HG",
    "AL",
    "GA",
    "CR",
    "V",
    "MO",
    "W",
}


class ProteinInspectionError(RuntimeError):
    """Raised when a protein structure cannot be inspected."""


def load_settings(settings_path: Path = SETTINGS_PATH) -> dict[str, Any]:
    """Load project settings from YAML."""

    if not settings_path.exists():
        raise FileNotFoundError(
            f"Settings file was not found: {settings_path}"
        )

    with settings_path.open("r", encoding="utf-8") as file:
        settings = yaml.safe_load(file)

    if not isinstance(settings, dict):
        raise ValueError("Settings file must contain a YAML dictionary.")

    return settings


def load_metadata(metadata_path: Path) -> dict[str, Any]:
    """Load RCSB metadata JSON."""

    if not metadata_path.exists():
        return {}

    try:
        with metadata_path.open("r", encoding="utf-8") as file:
            metadata = json.load(file)
    except json.JSONDecodeError as error:
        raise ProteinInspectionError(
            f"Invalid metadata JSON: {metadata_path}"
        ) from error

    if not isinstance(metadata, dict):
        return {}

    return metadata


def extract_entry_metadata(
    metadata: dict[str, Any],
) -> dict[str, Any]:
    """Extract useful experimental information from RCSB metadata."""

    title = (
        metadata.get("struct", {})
        .get("title")
    )

    experimental_methods = [
        item.get("method")
        for item in metadata.get("exptl", [])
        if item.get("method")
    ]

    resolution_values = (
        metadata.get("rcsb_entry_info", {})
        .get("resolution_combined", [])
    )

    deposit_date = (
        metadata.get("rcsb_accession_info", {})
        .get("deposit_date")
    )

    release_date = (
        metadata.get("rcsb_accession_info", {})
        .get("initial_release_date")
    )

    return {
        "title": title,
        "experimental_methods": experimental_methods,
        "resolution_angstrom": resolution_values,
        "deposit_date": deposit_date,
        "initial_release_date": release_date,
    }


def classify_hetero_residue(residue: Any) -> str:
    """Classify a non-protein residue."""

    residue_name = residue.get_resname().strip().upper()

    if residue_name in WATER_NAMES:
        return "water"

    if residue_name in COMMON_METALS:
        return "metal"

    atoms = list(residue.get_atoms())

    if (
        len(atoms) == 1
        and atoms[0].element.strip().upper() in COMMON_METALS
    ):
        return "metal"

    return "ligand_or_cofactor"


def inspect_structure(
    pdb_path: Path,
    metadata_path: Path | None = None,
) -> dict[str, Any]:
    """Inspect a PDB structure and return a structured report."""

    if not pdb_path.exists():
        raise FileNotFoundError(
            f"PDB structure was not found: {pdb_path}"
        )

    parser = PDBParser(QUIET=True)

    try:
        structure = parser.get_structure(
            pdb_path.stem,
            str(pdb_path),
        )
    except Exception as error:
        raise ProteinInspectionError(
            f"Could not parse PDB structure: {pdb_path}"
        ) from error

    models = list(structure.get_models())

    if not models:
        raise ProteinInspectionError(
            f"No structural models found in {pdb_path}"
        )

    first_model = models[0]

    total_atoms = 0
    protein_residue_count = 0
    water_count = 0

    chain_summaries: list[dict[str, Any]] = []
    ligands: list[dict[str, Any]] = []
    metals: list[dict[str, Any]] = []

    residue_name_counts: Counter[str] = Counter()
    ligand_name_counts: Counter[str] = Counter()
    metal_name_counts: Counter[str] = Counter()

    for chain in first_model:
        chain_atom_count = 0
        chain_protein_residues = 0
        chain_water_count = 0
        chain_hetero_count = 0

        for residue in chain:
            atoms = list(residue.get_atoms())
            atom_count = len(atoms)

            total_atoms += atom_count
            chain_atom_count += atom_count

            hetero_flag, residue_number, insertion_code = residue.id
            residue_name = residue.get_resname().strip().upper()

            if hetero_flag == " ":
                protein_residue_count += 1
                chain_protein_residues += 1
                residue_name_counts[residue_name] += 1
                continue

            classification = classify_hetero_residue(residue)

            residue_record = {
                "name": residue_name,
                "chain": chain.id,
                "residue_number": residue_number,
                "insertion_code": insertion_code.strip() or None,
                "atom_count": atom_count,
            }

            if classification == "water":
                water_count += 1
                chain_water_count += 1

            elif classification == "metal":
                metals.append(residue_record)
                metal_name_counts[residue_name] += 1
                chain_hetero_count += 1

            else:
                ligands.append(residue_record)
                ligand_name_counts[residue_name] += 1
                chain_hetero_count += 1

        chain_summaries.append(
            {
                "chain_id": chain.id,
                "protein_residues": chain_protein_residues,
                "atoms": chain_atom_count,
                "waters": chain_water_count,
                "hetero_residues": chain_hetero_count,
            }
        )

    metadata = {}

    if metadata_path is not None:
        metadata = extract_entry_metadata(
            load_metadata(metadata_path)
        )

    report = {
        "structure_id": pdb_path.stem.upper(),
        "source_file": str(pdb_path),
        "number_of_models": len(models),
        "number_of_chains": len(chain_summaries),
        "total_atoms": total_atoms,
        "protein_residues": protein_residue_count,
        "water_molecules": water_count,
        "ligand_instances": len(ligands),
        "metal_instances": len(metals),
        "chains": chain_summaries,
        "ligand_summary": dict(sorted(ligand_name_counts.items())),
        "metal_summary": dict(sorted(metal_name_counts.items())),
        "ligands": ligands,
        "metals": metals,
        "protein_residue_composition": dict(
            sorted(residue_name_counts.items())
        ),
        "experimental_metadata": metadata,
    }

    return report


def save_inspection_report(
    report: dict[str, Any],
) -> Path:
    """Save the inspection report as formatted JSON."""

    settings = load_settings()

    output_directory = (
        PROJECT_ROOT
        / settings["directories"]["protein_reports"]
    )

    output_directory.mkdir(parents=True, exist_ok=True)

    pdb_id = report["structure_id"]
    output_path = output_directory / f"{pdb_id}_inspection.json"

    with output_path.open("w", encoding="utf-8") as file:
        json.dump(report, file, indent=2)

    return output_path


def inspect_protein_by_id(pdb_id: str) -> tuple[dict[str, Any], Path]:
    """Inspect a locally downloaded protein using its PDB ID."""

    settings = load_settings()
    normalized_id = pdb_id.strip().upper()

    pdb_path = (
        PROJECT_ROOT
        / settings["directories"]["raw_proteins"]
        / f"{normalized_id}.pdb"
    )

    metadata_path = (
        PROJECT_ROOT
        / settings["directories"]["protein_metadata"]
        / f"{normalized_id}_metadata.json"
    )

    report = inspect_structure(
        pdb_path=pdb_path,
        metadata_path=metadata_path,
    )

    report_path = save_inspection_report(report)

    return report, report_path