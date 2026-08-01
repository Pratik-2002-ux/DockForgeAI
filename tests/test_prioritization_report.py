from pathlib import Path

from src.prioritization.prioritization_report import generate_prioritization_report


def test_generate_prioritization_report_writes_outputs(tmp_path: Path) -> None:
    report = generate_prioritization_report([{"compound_id": "A", "rank": 1, "consensus_score": 0.9}], tmp_path)

    assert (tmp_path / "consensus_ranking.csv").exists()
    assert (tmp_path / "diverse_top_hits.csv").exists()
    assert (tmp_path / "hit_prioritization_report.md").exists()
    assert (tmp_path / "hit_prioritization_report.json").exists()
    assert report["md_path"].endswith("hit_prioritization_report.md")
