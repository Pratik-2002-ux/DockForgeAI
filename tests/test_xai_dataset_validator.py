import pytest
from src.datasets.dataset_validator import DatasetValidationError, FINITE, MANDATORY, validate_records

def row(label="active",key="A",split="train"):
    value={column:"x" for column in MANDATORY}; value.update({column:1 for column in FINITE}); value.update({"compound_id":key,"canonical_isomeric_smiles":key,"inchikey":key,"binary_label":label,"dataset_split":split,"scaffold_id":key}); return value
def test_valid_records_pass() -> None:
    assert validate_records([row("active","A"),row("inactive_or_decoy","B")])["row_count"]==2
def test_duplicate_and_leakage_fail() -> None:
    duplicate=row("inactive_or_decoy","B"); duplicate["canonical_isomeric_smiles"]="A"; duplicate["scaffold_id"]="A"; duplicate["dataset_split"]="test"
    with pytest.raises(DatasetValidationError): validate_records([row("active","A"),duplicate])
