from pathlib import Path

from src.docking.report_generator import generate_reports


def test_generate_reports_writes_markdown_and_json(tmp_path: Path) -> None:
    analysis_csv = Path("data/processed/docking_tables/1HSG_MK1_redocking_analysis.csv")
    summary_path = Path("data/processed/docking_results/1HSG_MK1_redocking_summary.json")
    preparation_report_path = Path("data/processed/preparation_reports/1HSG_receptor_preparation.json")
    box_report_path = Path("data/processed/docking_boxes/1HSG_MK1_box.json")
    config_path = Path("config/docking/1HSG_MK1_redocking.txt")

    report_md_path = tmp_path / "report.md"
    report_json_path = tmp_path / "report.json"
    execution_metadata = {
        "executable": "vina",
        "vina_version": "1.2.7",
        "runtime_seconds": 42,
        "command": "vina --config config/docking/1HSG_MK1_redocking.txt",
    }

    report = generate_reports(
        analysis_csv=analysis_csv,
        summary_path=summary_path,
        execution_metadata=execution_metadata,
        report_md_path=report_md_path,
        report_json_path=report_json_path,
        preparation_report_path=preparation_report_path,
        box_report_path=box_report_path,
        config_path=config_path,
        threshold=2.0,
    )

    assert report_md_path.exists()
    assert report_json_path.exists()
    assert "# MK1 redocking report" in report_md_path.read_text(encoding="utf-8")
    assert report["summary"]["overall_pose_reproduction_success"] is True
