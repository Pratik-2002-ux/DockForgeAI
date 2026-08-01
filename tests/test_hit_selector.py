from pathlib import Path

from src.screening.hit_selector import select_hits


def test_select_hits_returns_top_candidates(tmp_path: Path) -> None:
    score_table = tmp_path / "scores.csv"
    score_table.write_text("compound_id,best_affinity,status,runtime_seconds,output_path\nA,-6.0,success,10,A\nB,-5.5,success,10,B\nC,-5.0,success,10,C\n", encoding="utf-8")
    output_path = tmp_path / "selected.csv"

    selected = select_hits(score_table, output_path, top_fraction=0.5, minimum_count=1, maximum_count=10)

    assert len(selected) == 2
    assert output_path.exists()
