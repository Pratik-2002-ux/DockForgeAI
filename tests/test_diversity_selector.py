from src.prioritization.diversity_selector import select_diverse_hits


def test_select_diverse_hits_returns_a_shortlist() -> None:
    rows = [{"compound_id": "A", "rank": 1, "consensus_score": 0.9}, {"compound_id": "B", "rank": 2, "consensus_score": 0.8}, {"compound_id": "C", "rank": 3, "consensus_score": 0.7}]

    selected = select_diverse_hits(rows, top_n=2)

    assert len(selected) == 2
