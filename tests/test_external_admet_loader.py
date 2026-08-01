from pathlib import Path

from src.admet.external_admet_loader import load_external_admet


def test_load_external_admet_reads_csv(tmp_path: Path) -> None:
    table = tmp_path / "admet.csv"
    table.write_text("compound_id,source_tool\nA,testtool\n", encoding="utf-8")

    rows = load_external_admet(table)

    assert rows[0]["compound_id"] == "A"
