"""Deterministic 3D preparation of the committed Top-30 review pool."""
from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from src.docking.phase4_common import atomic_json, provenance, sha256, write_checksums, write_reports

REQUIRED = {"compound_id", "canonical_isomeric_smiles", "supplier_id", "catalogue_id", "availability", "predicted_pki", "applicability_domain_status", "pains_flag", "brenk_flag"}


class LigandPreparationError(RuntimeError):
    """Raised when a ligand cannot be prepared without inventing data."""


def _run(command: list[str], label: str) -> None:
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    if completed.returncode:
        raise LigandPreparationError(f"{label} failed: {completed.stderr.strip() or completed.stdout.strip()}")


def _prepare_one(row: pd.Series, output_dir: Path, seed: int) -> dict[str, Any]:
    compound_id = str(row["compound_id"])
    molecule = Chem.MolFromSmiles(str(row["canonical_isomeric_smiles"]))
    if molecule is None:
        raise LigandPreparationError(f"Invalid SMILES: {compound_id}")
    canonical = Chem.MolToSmiles(molecule, isomericSmiles=True, canonical=True)
    molecule = Chem.AddHs(molecule)
    parameters = AllChem.ETKDGv3(); parameters.randomSeed = seed; parameters.numThreads = 1
    if AllChem.EmbedMolecule(molecule, parameters) != 0:
        raise LigandPreparationError(f"3D embedding failed: {compound_id}")
    forcefield = "MMFF94"
    properties = AllChem.MMFFGetMoleculeProperties(molecule, mmffVariant="MMFF94")
    if properties is not None:
        field = AllChem.MMFFGetMoleculeForceField(molecule, properties)
    else:
        forcefield = "UFF"; field = AllChem.UFFGetMoleculeForceField(molecule)
    if field is None:
        raise LigandPreparationError(f"No supported force field: {compound_id}")
    energy_before = float(field.CalcEnergy()); convergence = int(field.Minimize(maxIts=1000)); energy_after = float(field.CalcEnergy())
    molecule.SetProp("compound_id", compound_id); molecule.SetProp("canonical_isomeric_smiles", canonical)
    sdf = output_dir / "sdf" / f"{compound_id}.sdf"; mol2 = output_dir / "mol2" / f"{compound_id}.mol2"; pdbqt = output_dir / "pdbqt" / f"{compound_id}.pdbqt"
    sdf.parent.mkdir(parents=True, exist_ok=True); mol2.parent.mkdir(parents=True, exist_ok=True); pdbqt.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(sdf)); writer.write(molecule); writer.close()
    obabel = shutil.which("obabel")
    meeko = shutil.which("mk_prepare_ligand.py")
    if not obabel:
        raise LigandPreparationError("Open Babel executable 'obabel' is unavailable")
    if not meeko:
        raise LigandPreparationError("Meeko executable 'mk_prepare_ligand.py' is unavailable")
    _run([obabel, str(sdf), "-O", str(mol2)], f"MOL2 conversion for {compound_id}")
    _run([meeko, "-i", str(sdf), "-o", str(pdbqt)], f"PDBQT preparation for {compound_id}")
    if not all(path.is_file() and path.stat().st_size for path in (sdf, mol2, pdbqt)):
        raise LigandPreparationError(f"Incomplete output for {compound_id}")
    return {"compound_id": compound_id, "canonical_isomeric_smiles": canonical, "forcefield": forcefield, "minimization_convergence_code": convergence, "energy_before": energy_before, "energy_after": energy_after, "protonation_method": "RDKit explicit hydrogens; Meeko Gasteiger preparation", "sdf_path": str(sdf), "mol2_path": str(mol2), "pdbqt_path": str(pdbqt), "sdf_checksum": sha256(sdf), "mol2_checksum": sha256(mol2), "pdbqt_checksum": sha256(pdbqt)}


def prepare_ligands(input_path: Path, output_dir: Path, root: Path, overwrite: bool = False, seed: int = 2026) -> dict[str, Any]:
    """Prepare every ligand in a committed review-pool CSV."""

    if output_dir.exists() and any(output_dir.iterdir()) and not overwrite:
        raise FileExistsError("Ligand outputs exist; use --overwrite")
    if overwrite and output_dir.exists():
        for child in output_dir.iterdir():
            if child.is_dir(): shutil.rmtree(child)
            else: child.unlink()
    output_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(input_path)
    missing = sorted(REQUIRED - set(frame.columns))
    if missing:
        raise ValueError(f"Review pool missing required columns: {missing}")
    screening_checks = json.loads((input_path.parent / "screening_checksums.json").read_text())
    if input_path.name in screening_checks and sha256(input_path) != screening_checks[input_path.name]:
        raise ValueError("Top-30 review-pool checksum mismatch")
    records: list[dict[str, Any]] = []
    for index, (_, row) in enumerate(frame.iterrows(), 1):
        print(f"Preparing ligand {index}/{len(frame)}: {row['compound_id']}", flush=True)
        records.append(_prepare_one(row, output_dir, seed))
    prepared_table = pd.DataFrame(records)
    prepared_table.to_csv(output_dir / "prepared_ligands.csv", index=False)
    aggregate = Chem.SDWriter(str(output_dir / "prepared_ligands.sdf"))
    for record in records:
        supplier = Chem.SDMolSupplier(record["sdf_path"], removeHs=False)[0]
        if supplier is not None: aggregate.write(supplier)
    aggregate.close()
    displayed_input = str(input_path.relative_to(root)) if input_path.is_relative_to(root) else str(input_path)
    manifest = {"schema_version": 1, "module": "4.2_ligand_preparation", "input_path": displayed_input, "input_checksum": sha256(input_path), "input_count": len(frame), "prepared_count": len(records), "failed_count": 0, "conformer_method": "RDKit ETKDGv3", "minimization": "MMFF94 with UFF fallback", "protonation": "explicit hydrogens plus Meeko ligand preparation; no experimental microspecies assignment", **provenance(root, seed)}
    report_values = {"Input ligands": len(frame), "Prepared ligands": len(records), "3D method": manifest["conformer_method"], "Minimization": manifest["minimization"], "Random seed": seed}
    write_reports("DockForgeAI Ligand Preparation", report_values, output_dir / "ligand_report.html", output_dir / "ligand_report.pdf")
    data_outputs = [output_dir / "prepared_ligands.csv", output_dir / "prepared_ligands.sdf", output_dir / "ligand_report.html", output_dir / "ligand_report.pdf", *[Path(record[key]) for record in records for key in ("sdf_path", "mol2_path", "pdbqt_path")]]
    manifest["output_checksums"] = {str(path.relative_to(output_dir)): sha256(path) for path in data_outputs}
    atomic_json(manifest, output_dir / "ligand_manifest.json")
    write_checksums([*data_outputs, output_dir / "ligand_manifest.json"], output_dir / "ligand_checksums.json")
    return manifest
