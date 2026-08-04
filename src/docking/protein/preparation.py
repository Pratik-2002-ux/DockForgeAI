"""Deterministic, provenance-rich receptor preparation for Phase 4."""
from __future__ import annotations

import shutil
import urllib.request
from pathlib import Path
from typing import Any

from Bio import PDB

from src.docking.phase4_common import atomic_json, provenance, sha256, write_checksums, write_reports

WATERS = {"HOH", "WAT", "DOD"}


def download_pdb(pdb_id: str, destination: Path, allow_network: bool = False) -> Path:
    """Return a cached PDB, downloading from RCSB only when explicitly allowed."""

    if destination.is_file() and destination.stat().st_size:
        return destination
    if not allow_network:
        raise FileNotFoundError(f"Cached PDB unavailable and network download disabled: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb_id}.pdb", destination)
    return destination


def validate_structure(path: Path) -> dict[str, Any]:
    """Validate chains and catalytic Asp25 residues in a receptor PDB."""

    structure = PDB.PDBParser(QUIET=True).get_structure("receptor", str(path))
    model = next(structure.get_models())
    chains = [chain.id for chain in model]
    catalytic = [f"{chain.id}:25" for chain in model for residue in chain if residue.id[0] == " " and residue.id[1] == 25 and residue.get_resname().strip() == "ASP"]
    if not {"A", "B"}.issubset(chains) or not {"A:25", "B:25"}.issubset(catalytic):
        raise ValueError("Receptor must contain chains A/B and catalytic ASP25 in both chains")
    atoms = sum(1 for _ in structure.get_atoms())
    residues = sum(1 for _ in structure.get_residues())
    return {"chains": chains, "catalytic_residues": catalytic, "atom_count": atoms, "residue_count": residues}


class _ResidueFilter(PDB.Select):
    def __init__(self, remove_water: bool, remove_hetero: bool) -> None:
        self.remove_water = remove_water
        self.remove_hetero = remove_hetero

    def accept_residue(self, residue: PDB.Residue.Residue) -> bool:
        name = residue.get_resname().strip().upper()
        if self.remove_water and name in WATERS:
            return False
        if self.remove_hetero and residue.id[0] != " " and name not in WATERS:
            return False
        return True


def _filter_pdb(source: Path, destination: Path, remove_water: bool, remove_hetero: bool) -> Path:
    structure = PDB.PDBParser(QUIET=True).get_structure("receptor", str(source))
    destination.parent.mkdir(parents=True, exist_ok=True)
    writer = PDB.PDBIO(); writer.set_structure(structure); writer.save(str(destination), _ResidueFilter(remove_water, remove_hetero))
    return destination


def remove_waters(source: Path, destination: Path) -> Path:
    """Remove crystallographic waters."""

    return _filter_pdb(source, destination, True, False)


def remove_ligands(source: Path, destination: Path) -> Path:
    """Remove non-water hetero residues."""

    return _filter_pdb(source, destination, False, True)


def repair_missing_atoms(source: Path, destination: Path) -> dict[str, Any]:
    """Preserve an already curated receptor; do not invent coordinates."""

    shutil.copy2(source, destination)
    return {"path": str(destination), "method": "curated_cached_structure_no_coordinate_invention", "atoms_added": 0}


def assign_protonation(prepared_pdb: Path, prepared_pdbqt: Path, cached_pdb: Path, cached_pdbqt: Path) -> dict[str, Any]:
    """Use the previously validated Meeko-prepared receptor artifacts."""

    if not cached_pdb.is_file() or not cached_pdbqt.is_file():
        raise FileNotFoundError("Validated cached Meeko receptor PDB/PDBQT artifacts are required")
    shutil.copy2(cached_pdb, prepared_pdb); shutil.copy2(cached_pdbqt, prepared_pdbqt)
    return {"method": "validated_cached_meeko_default_templates", "assumption": "standard Asp carboxylate templates; protonation was not experimentally assigned"}


def save_manifest(manifest: dict[str, Any], path: Path) -> None:
    """Write the protein manifest atomically."""

    atomic_json(manifest, path)


def prepare_receptor(root: Path, output_dir: Path, pdb_id: str = "1HSG", overwrite: bool = False) -> dict[str, Any]:
    """Prepare the cached HIV-1 protease receptor and all required reports."""

    names = ["receptor_clean.pdb", "receptor_prepared.pdb", "receptor_prepared.pdbqt", "protein_manifest.json", "protein_checksums.json", "protein_report.html", "protein_report.pdf"]
    if not overwrite and any((output_dir / name).exists() for name in names):
        raise FileExistsError("Protein outputs exist; use --overwrite")
    output_dir.mkdir(parents=True, exist_ok=True)
    source = download_pdb(pdb_id, root / "data/raw/proteins/1HSG.pdb")
    waterless = output_dir / "receptor_without_waters.pdb"
    clean = output_dir / "receptor_clean.pdb"
    remove_waters(source, waterless); remove_ligands(waterless, clean)
    repair = repair_missing_atoms(clean, output_dir / "receptor_repaired.pdb")
    protonation = assign_protonation(output_dir / "receptor_prepared.pdb", output_dir / "receptor_prepared.pdbqt", root / "data/processed/protonated_receptors/1HSG_chains-A-B_receptor.pdb", root / "data/processed/pdbqt/receptors/1HSG_chains-A-B_receptor.pdbqt")
    validation = validate_structure(output_dir / "receptor_prepared.pdb")
    manifest = {"schema_version": 1, "module": "4.1_protein_preparation", "pdb_id": pdb_id, "source_path": str(source.relative_to(root)), "source_checksum": sha256(source), "validation": validation, "repair": repair, "protonation": protonation, **provenance(root)}
    report_values = {"PDB ID": pdb_id, "Chains": ", ".join(validation["chains"]), "Catalytic residues": ", ".join(validation["catalytic_residues"]), "Atoms": validation["atom_count"], "Protonation": protonation["method"]}
    write_reports("DockForgeAI Protein Preparation", report_values, output_dir / "protein_report.html", output_dir / "protein_report.pdf")
    save_manifest(manifest, output_dir / "protein_manifest.json")
    data_artifacts = [output_dir / name for name in ["receptor_clean.pdb", "receptor_repaired.pdb", "receptor_prepared.pdb", "receptor_prepared.pdbqt", "protein_report.html", "protein_report.pdf"]]
    manifest["output_checksums"] = {path.name: sha256(path) for path in data_artifacts}
    save_manifest(manifest, output_dir / "protein_manifest.json")
    write_checksums([*data_artifacts, output_dir / "protein_manifest.json"], output_dir / "protein_checksums.json")
    return manifest
