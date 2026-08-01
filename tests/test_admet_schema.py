from pathlib import Path

import pytest

from src.admet.admet_schema import ADMETSchemaError, validate_admet_table


def test_validate_admet_table_accepts_valid_rows(tmp_path: Path) -> None:
    table = tmp_path / "admet.csv"
    table.write_text("compound_id,source_tool,ames_probability\nA,tool,0.2\n", encoding="utf-8")

    result = validate_admet_table(table)

    assert result["row_count"] == 1


def test_validate_admet_table_rejects_invalid_probability(tmp_path: Path) -> None:
    table = tmp_path / "admet.csv"
    table.write_text("compound_id,source_tool,ames_probability\nA,tool,1.2\n", encoding="utf-8")

    with pytest.raises(ADMETSchemaError):
        validate_admet_table(table)
