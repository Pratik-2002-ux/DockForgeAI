"""Tests for Phase 4 protein preparation."""
from pathlib import Path

import pytest

from src.docking.protein import download_pdb, prepare_receptor, validate_structure

ROOT = Path(__file__).resolve().parents[1]


def test_validate_structure() -> None:
    result = validate_structure(ROOT / "data/processed/protonated_receptors/1HSG_chains-A-B_receptor.pdb")
    assert {"A:25", "B:25"}.issubset(result["catalytic_residues"])


def test_download_requires_explicit_network(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="network download disabled"):
        download_pdb("1HSG", tmp_path / "missing.pdb")


def test_prepare_protein_outputs(tmp_path: Path) -> None:
    manifest = prepare_receptor(ROOT, tmp_path / "protein")
    assert manifest["pdb_id"] == "1HSG"
    for name in ["protein_manifest.json", "protein_checksums.json", "protein_report.html", "protein_report.pdf", "receptor_prepared.pdbqt"]:
        assert (tmp_path / "protein" / name).stat().st_size > 0


def test_prepare_protein_overwrite_guard(tmp_path: Path) -> None:
    output = tmp_path / "protein"; prepare_receptor(ROOT, output)
    with pytest.raises(FileExistsError):
        prepare_receptor(ROOT, output)
