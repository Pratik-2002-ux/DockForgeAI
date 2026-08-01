from src.ligands.druglikeness import calculate_descriptors


def test_calculate_descriptors_returns_expected_fields() -> None:
    records = [{"compound_id": "A", "smiles": "CCO"}]

    descriptors = calculate_descriptors(records)

    assert descriptors[0]["compound_id"] == "A"
    assert descriptors[0]["molecular_weight"] > 0.0
    assert descriptors[0]["qed"] >= 0.0
