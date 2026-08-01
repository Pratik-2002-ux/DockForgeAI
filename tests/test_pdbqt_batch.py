from pathlib import Path

from src.ligands.pdbqt_batch import prepare_pdbqt


def test_prepare_pdbqt_skips_existing_outputs(tmp_path: Path) -> None:
    output_dir = tmp_path / "pdbqt"
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "A.pdbqt").write_text("ROOT\nENDROOT\nTORSDOF\n", encoding="utf-8")

    result = prepare_pdbqt([
        {"compound_id": "A", "output_path": str(tmp_path / "A.sdf"), "status": "ready"}
    ], output_dir=output_dir)

    assert result["prepared"][0]["status"] == "skipped"
