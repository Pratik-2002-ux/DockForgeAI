from src.prioritization.feature_normalizer import normalize_features


def test_normalize_features_handles_beneficial_and_risk_metrics() -> None:
    rows = [{"compound_id": "A", "qed": 0.5, "best_affinity": -6.0}, {"compound_id": "B", "qed": 0.8, "best_affinity": -4.0}]

    normalized = normalize_features(rows, beneficial={"qed"}, risk={"best_affinity"})

    assert normalized[0]["norm_qed"] >= 0.0
    assert normalized[0]["norm_best_affinity"] >= 0.0
