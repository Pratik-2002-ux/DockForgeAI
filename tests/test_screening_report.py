from pathlib import Path

from src.screening.screening_report import generate_screening_report


def test_generate_screening_report_writes_outputs(tmp_path: Path) -> None:
    report = generate_screening_report({"library_size": 3, "accepted_count": 2, "rejected_count": 1}, tmp_path)

    assert (tmp_path / "screening_report.md").exists()
    assert (tmp_path / "screening_report.json").exists()
    assert report["md_path"].endswith("screening_report.md")
