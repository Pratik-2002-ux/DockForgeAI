"""Tests for deterministic Phase 4 ligand preparation."""
from pathlib import Path

import pandas as pd
import pytest

from src.docking.ligands import prepare_ligands

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/processed/prospective_studies/HIVPR_2026_001/ml_screening/top_30_diverse_review_pool.csv"


def test_prepare_ligands_outputs(tmp_path: Path) -> None:
    source = tmp_path / "top_30_diverse_review_pool.csv"
    pd.read_csv(SOURCE).head(1).to_csv(source, index=False)
    (tmp_path / "screening_checksums.json").write_text("{}")
    manifest = prepare_ligands(source, tmp_path / "ligands", ROOT)
    assert manifest["prepared_count"] == 1
    for suffix in ("sdf", "mol2", "pdbqt"):
        assert len(list((tmp_path / "ligands" / suffix).glob(f"*.{suffix}"))) == 1


def test_missing_columns_fail(tmp_path: Path) -> None:
    source = tmp_path / "top_30_diverse_review_pool.csv"; pd.DataFrame({"compound_id": ["x"]}).to_csv(source, index=False)
    (tmp_path / "screening_checksums.json").write_text("{}")
    with pytest.raises(ValueError, match="missing required columns"):
        prepare_ligands(source, tmp_path / "ligands", ROOT)


def test_ligand_overwrite_guard(tmp_path: Path) -> None:
    source = tmp_path / "top_30_diverse_review_pool.csv"; pd.read_csv(SOURCE).head(1).to_csv(source, index=False)
    (tmp_path / "screening_checksums.json").write_text("{}")
    output = tmp_path / "ligands"; prepare_ligands(source, output, ROOT)
    with pytest.raises(FileExistsError): prepare_ligands(source, output, ROOT)
