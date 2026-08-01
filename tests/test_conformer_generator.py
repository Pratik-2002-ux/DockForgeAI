from pathlib import Path

from src.ligands.conformer_generator import generate_conformers


def test_generate_conformers_writes_sdf_files(tmp_path: Path) -> None:
    records = [{"compound_id": "A", "smiles": "CCO"}]

    generated = generate_conformers(records, tmp_path)

    assert len(generated) == 1
    assert (tmp_path / "A.sdf").exists()
