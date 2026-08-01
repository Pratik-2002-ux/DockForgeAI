from pathlib import Path

from src.ligands.library_loader import load_library


def test_load_library_from_csv(tmp_path: Path) -> None:
    library_path = tmp_path / "library.csv"
    library_path.write_text("compound_id,smiles\nA,CCO\nB,CCN\n", encoding="utf-8")

    summary = load_library(library_path, max_compounds=5)

    assert summary["record_count"] == 2
    assert summary["records"][0]["compound_id"] == "A"
