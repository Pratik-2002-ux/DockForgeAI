from pathlib import Path

from src.screening.score_parser import parse_scores


def test_parse_scores_writes_csv(tmp_path: Path) -> None:
    output_path = tmp_path / "scores.csv"
    rows = parse_scores([
        {"compound_id": "A", "best_affinity": -6.0, "status": "success", "output_path": "A.pdbqt"}
    ], output_path)

    assert rows[0]["compound_id"] == "A"
    assert output_path.exists()
