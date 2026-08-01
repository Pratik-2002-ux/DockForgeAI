from pathlib import Path

import pytest

from src.analysis.complex_builder import ComplexBuilderError, build_complex


def test_build_complex_writes_receptor_and_ligand_records(tmp_path: Path) -> None:
    receptor = Path("data/processed/protonated_receptors/1HSG_chains-A-B_receptor.pdb")
    ligand = Path("data/processed/crystallographic_ligands/1HSG_MK1_crystallographic_corrected.sdf")
    output = tmp_path / "complex.pdb"

    result = build_complex(receptor_path=receptor, ligand_path=ligand, output_path=output, label="crystal")

    text = output.read_text(encoding="utf-8")
    assert result["contains_receptor"] is True
    assert result["contains_ligand"] is True
    assert "ATOM" in text
    assert "TER" in text


def test_build_complex_rejects_missing_input(tmp_path: Path) -> None:
    with pytest.raises(ComplexBuilderError):
        build_complex(receptor_path=tmp_path / "missing.pdb", ligand_path=tmp_path / "missing.sdf", output_path=tmp_path / "complex.pdb", label="crystal")
