from src.ligands.standardizer import standardize_ligands


def test_standardize_ligands_deduplicates_and_keeps_valid_records() -> None:
    records = [
        {"compound_id": "A", "smiles": "CCO", "source": "demo"},
        {"compound_id": "B", "smiles": "CCO", "source": "demo"},
        {"compound_id": "C", "smiles": "CCN", "source": "demo"},
    ]

    result = standardize_ligands(records)

    assert len(result["standardized"]) == 2
    assert len(result["rejected"]) == 1
