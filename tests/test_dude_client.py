from pathlib import Path
import pytest
from src.datasets.dude_client import DUDEClient, DUDEClientError

def test_parse_inspected_ism(tmp_path: Path) -> None:
    path=tmp_path/"decoys_final.ism"; path.write_text("CCO ZINC1\nCCC ZINC2 extra\n")
    assert [r["source_molecule_id"] for r in DUDEClient.parse_ism(path)] == ["ZINC1","ZINC2"]
def test_malformed_ism_rejected(tmp_path: Path) -> None:
    path=tmp_path/"bad.ism"; path.write_text("CCO\n")
    with pytest.raises(DUDEClientError): DUDEClient.parse_ism(path)
