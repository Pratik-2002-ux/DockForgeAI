from src.prioritization.consensus_ranker import compute_consensus_ranking


def test_compute_consensus_ranking_ranks_records() -> None:
    rows = [{"compound_id": "A", "norm_best_affinity": 0.7, "norm_interaction": 0.8, "norm_qed": 0.6, "norm_alerts": 0.9, "norm_ligand_efficiency": 0.5}, {"compound_id": "B", "norm_best_affinity": 0.9, "norm_interaction": 0.7, "norm_qed": 0.8, "norm_alerts": 0.8, "norm_ligand_efficiency": 0.7}]

    ranked = compute_consensus_ranking(rows, {"docking": 0.30, "interaction": 0.25, "druglikeness": 0.20, "structural_alerts": 0.10, "ligand_efficiency": 0.15})

    assert ranked[0]["rank"] == 1
    assert ranked[0]["compound_id"] in {"A", "B"}
