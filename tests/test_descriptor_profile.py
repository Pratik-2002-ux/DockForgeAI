from src.admet.descriptor_profile import build_descriptor_profile


def test_build_descriptor_profile_returns_profiles() -> None:
    records = [{"compound_id": "A", "smiles": "CCO", "best_affinity": -6.0}]

    profiles = build_descriptor_profile(records)

    assert profiles[0]["compound_id"] == "A"
    assert profiles[0]["molecular_weight"] > 0.0
