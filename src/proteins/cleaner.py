"""Protein cleaning and co-crystallized ligand extraction.

This module prepares an intermediate receptor structure by:

1. Preserving selected protein chains.
2. Removing crystallographic water molecules.
3. Extracting selected ligands into separate PDB files.
4. Removing remaining non-protein hetero-residues unless retained.
5. Writing a preparation report for reproducibility.

This stage does not yet add hydrogens, assign protonation states,
calculate charges or generate PDBQT files.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import yaml
from Bio.PDB import PDBIO, PDBParser, Select


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SETTINGS_PATH = PROJECT_ROOT / "config" / "settings.yaml"

WATER_NAMES = {
    "HOH",
    "WAT",
    "DOD",
}


class ProteinCleaningError(RuntimeError):
    """Raised when protein cleaning cannot be completed."""


def load_settings(
    settings_path: Path = SETTINGS_PATH,
) -> dict[str, Any]:
    """Load YAML settings."""

    if not settings_path.exists():
        raise FileNotFoundError(
            f"Settings file was not found: {settings_path}"
        )

    with settings_path.open("r", encoding="utf-8") as file:
        settings = yaml.safe_load(file)

    if not isinstance(settings, dict):
        raise ValueError(
            "Settings file must contain a YAML dictionary."
        )

    return settings


def normalize_names(names: list[str] | None) -> set[str]:
    """Normalize residue or ligand names."""

    if not names:
        return set()

    return {
        name.strip().upper()
        for name in names
        if name.strip()
    }


class ReceptorSelect(Select):
    """Select protein residues and explicitly retained hetero-residues."""

    def __init__(
        self,
        selected_chains: set[str],
        retained_hetero_names: set[str],
    ) -> None:
        self.selected_chains = selected_chains
        self.retained_hetero_names = retained_hetero_names

    def accept_chain(self, chain: Any) -> int:
        """Keep only selected chains."""

        return int(chain.id in self.selected_chains)

    def accept_residue(self, residue: Any) -> int:
        """Keep protein residues and selected hetero-residues."""

        hetero_flag = residue.id[0]
        residue_name = residue.get_resname().strip().upper()

        if hetero_flag == " ":
            return 1

        if residue_name in WATER_NAMES:
            return 0

        return int(residue_name in self.retained_hetero_names)


class LigandSelect(Select):
    """Select a specific ligand from a specific chain and residue number."""

    def __init__(
        self,
        ligand_name: str,
        chain_id: str,
        residue_number: int,
    ) -> None:
        self.ligand_name = ligand_name.upper()
        self.chain_id = chain_id
        self.residue_number = residue_number

    def accept_chain(self, chain: Any) -> int:
        return int(chain.id == self.chain_id)

    def accept_residue(self, residue: Any) -> int:
        residue_name = residue.get_resname().strip().upper()
        residue_number = residue.id[1]

        return int(
            residue_name == self.ligand_name
            and residue_number == self.residue_number
        )


def inspect_cleaning_candidates(
    structure: Any,
    selected_chains: set[str],
) -> dict[str, Any]:
    """Inspect selected chains before cleaning."""

    protein_residues = 0
    waters = 0
    hetero_residues: list[dict[str, Any]] = []

    model = next(structure.get_models())

    available_chains = {
        chain.id
        for chain in model
    }

    missing_chains = selected_chains - available_chains

    if missing_chains:
        raise ProteinCleaningError(
            "Requested chain IDs were not found: "
            + ", ".join(sorted(missing_chains))
        )

    for chain in model:
        if chain.id not in selected_chains:
            continue

        for residue in chain:
            hetero_flag, residue_number, insertion_code = residue.id
            residue_name = residue.get_resname().strip().upper()

            if hetero_flag == " ":
                protein_residues += 1
                continue

            if residue_name in WATER_NAMES:
                waters += 1
                continue

            hetero_residues.append(
                {
                    "name": residue_name,
                    "chain": chain.id,
                    "residue_number": residue_number,
                    "insertion_code": (
                        insertion_code.strip() or None
                    ),
                    "atom_count": len(list(residue.get_atoms())),
                }
            )

    return {
        "selected_chains": sorted(selected_chains),
        "protein_residues": protein_residues,
        "waters": waters,
        "hetero_residues": hetero_residues,
    }


def extract_selected_ligands(
    structure: Any,
    ligand_names: set[str],
    selected_chains: set[str],
    output_directory: Path,
    pdb_id: str,
) -> list[dict[str, Any]]:
    """Extract matching ligands into individual PDB files."""

    extracted: list[dict[str, Any]] = []
    model = next(structure.get_models())

    output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    for chain in model:
        if chain.id not in selected_chains:
            continue

        for residue in chain:
            hetero_flag, residue_number, insertion_code = residue.id
            residue_name = residue.get_resname().strip().upper()

            if hetero_flag == " ":
                continue

            if residue_name in WATER_NAMES:
                continue

            if residue_name not in ligand_names:
                continue

            output_name = (
                f"{pdb_id}_{residue_name}_"
                f"chain-{chain.id}_res-{residue_number}.pdb"
            )

            output_path = output_directory / output_name

            io = PDBIO()
            io.set_structure(structure)
            io.save(
                str(output_path),
                LigandSelect(
                    ligand_name=residue_name,
                    chain_id=chain.id,
                    residue_number=residue_number,
                ),
            )

            extracted.append(
                {
                    "name": residue_name,
                    "chain": chain.id,
                    "residue_number": residue_number,
                    "insertion_code": (
                        insertion_code.strip() or None
                    ),
                    "output_file": str(output_path),
                    "atom_count": len(list(residue.get_atoms())),
                }
            )

    return extracted


def clean_protein(
    pdb_id: str,
    selected_chains: list[str],
    extract_ligands: list[str] | None = None,
    retain_hetero: list[str] | None = None,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Clean a downloaded protein structure.

    Parameters
    ----------
    pdb_id:
        Four-character PDB identifier.
    selected_chains:
        Protein chains to retain.
    extract_ligands:
        Ligand residue names to save separately and remove from receptor.
    retain_hetero:
        Hetero-residue names to preserve in the receptor, for example
        essential cofactors or catalytic metal-containing groups.
    overwrite:
        Replace existing output files.
    """

    settings = load_settings()
    normalized_id = pdb_id.strip().upper()

    chain_ids = {
        chain.strip()
        for chain in selected_chains
        if chain.strip()
    }

    if not chain_ids:
        raise ValueError(
            "At least one protein chain must be selected."
        )

    ligand_names = normalize_names(extract_ligands)
    retained_names = normalize_names(retain_hetero)

    overlap = ligand_names & retained_names

    if overlap:
        raise ValueError(
            "A residue cannot be both extracted and retained: "
            + ", ".join(sorted(overlap))
        )

    input_path = (
        PROJECT_ROOT
        / settings["directories"]["raw_proteins"]
        / f"{normalized_id}.pdb"
    )

    if not input_path.exists():
        raise FileNotFoundError(
            f"Downloaded PDB file was not found: {input_path}"
        )

    receptor_directory = (
        PROJECT_ROOT
        / settings["directories"]["prepared_proteins"]
    )

    ligand_directory = (
        PROJECT_ROOT
        / settings["directories"]["extracted_ligands"]
    )

    report_directory = (
        PROJECT_ROOT
        / settings["directories"]["preparation_reports"]
    )

    receptor_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    ligand_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    report_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    chain_label = "-".join(sorted(chain_ids))

    receptor_path = (
        receptor_directory
        / f"{normalized_id}_chains-{chain_label}_clean.pdb"
    )

    report_path = (
        report_directory
        / f"{normalized_id}_cleaning_report.json"
    )

    if receptor_path.exists() and not overwrite:
        raise FileExistsError(
            f"Clean receptor already exists: {receptor_path}. "
            "Use --overwrite to replace it."
        )

    parser = PDBParser(QUIET=True)

    try:
        structure = parser.get_structure(
            normalized_id,
            str(input_path),
        )
    except Exception as error:
        raise ProteinCleaningError(
            f"Could not parse protein structure: {input_path}"
        ) from error

    candidate_summary = inspect_cleaning_candidates(
        structure=structure,
        selected_chains=chain_ids,
    )

    extracted = extract_selected_ligands(
        structure=structure,
        ligand_names=ligand_names,
        selected_chains=chain_ids,
        output_directory=ligand_directory,
        pdb_id=normalized_id,
    )

    found_ligand_names = {
        item["name"]
        for item in extracted
    }

    missing_ligands = ligand_names - found_ligand_names

    if missing_ligands:
        raise ProteinCleaningError(
            "Requested ligand residue names were not found in "
            "the selected chains: "
            + ", ".join(sorted(missing_ligands))
        )

    io = PDBIO()
    io.set_structure(structure)
    io.save(
        str(receptor_path),
        ReceptorSelect(
            selected_chains=chain_ids,
            retained_hetero_names=retained_names,
        ),
    )

    if not receptor_path.exists() or receptor_path.stat().st_size == 0:
        raise ProteinCleaningError(
            "The cleaned receptor file was not created correctly."
        )

    removed_hetero = [
        item
        for item in candidate_summary["hetero_residues"]
        if item["name"] not in retained_names
    ]

    retained_hetero_records = [
        item
        for item in candidate_summary["hetero_residues"]
        if item["name"] in retained_names
    ]

    report = {
        "pdb_id": normalized_id,
        "input_file": str(input_path),
        "clean_receptor_file": str(receptor_path),
        "selected_chains": sorted(chain_ids),
        "protein_residues_retained": (
            candidate_summary["protein_residues"]
        ),
        "waters_removed": candidate_summary["waters"],
        "ligands_requested_for_extraction": sorted(ligand_names),
        "extracted_ligands": extracted,
        "retained_hetero_residues": retained_hetero_records,
        "removed_hetero_residues": removed_hetero,
        "important_note": (
            "This is an intermediate cleaned PDB. Hydrogens, "
            "protonation states, partial charges and PDBQT preparation "
            "have not yet been performed."
        ),
    }

    with report_path.open("w", encoding="utf-8") as file:
        json.dump(
            report,
            file,
            indent=2,
        )

    report["report_file"] = str(report_path)

    return report