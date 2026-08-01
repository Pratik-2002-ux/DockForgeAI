from pathlib import Path

from src.docking.redocking_validator import validate_redocking


def test_validate_redocking_summarizes_success(tmp_path: Path) -> None:
    analysis_csv = Path("data/processed/docking_tables/1HSG_MK1_redocking_analysis.csv")
    summary_path = tmp_path / "summary.json"

    summary = validate_redocking(analysis_csv=analysis_csv, summary_path=summary_path, threshold=2.0)

    assert summary["overall_pose_reproduction_success"] is True
    assert summary["top_ranked_pose_success"] is True
    assert summary["passing_pose_percentage"] == 10.0
    assert summary["first_passing_pose_rank"] == 1
    assert summary["min_rmsd_pose_rank"] == 1
    assert summary["min_rmsd_angstrom"] == 0.518
    assert summary["passing_pose_count"] == 1
    assert summary["interpretation"] == "Validated strongly"
