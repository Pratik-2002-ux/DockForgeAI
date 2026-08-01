from pathlib import Path
from src.reports.report_data import collect_report_data
from src.reports.report_figures import generate_report_figures

def test_report_figures_created_and_nonempty(tmp_path: Path) -> None:
    outputs = generate_report_figures(collect_report_data(), tmp_path)
    assert len(outputs) >= 7 and all(path.stat().st_size > 0 for path in outputs)

def test_missing_figure_sections_skip_cleanly(tmp_path: Path) -> None:
    assert generate_report_figures({"redocking": {}, "interactions": {}, "screening": {}, "prioritization": {}}, tmp_path) == []
