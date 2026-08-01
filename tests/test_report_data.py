from pathlib import Path
import yaml
from src.reports.report_data import collect_report_data

def test_consolidated_report_data_has_required_sections() -> None:
    data = collect_report_data()
    assert {"project", "software", "protein", "redocking", "interactions", "screening", "prioritization", "limitations", "files"} <= data.keys()
    assert data["redocking"]["rmsd_angstrom"] == 0.518

def test_missing_optional_outputs_do_not_crash(tmp_path: Path) -> None:
    (tmp_path / "config").mkdir(); path = tmp_path / "config/settings.yaml"
    path.write_text(yaml.safe_dump({"project": {"name": "test", "version": "1"}, "directories": {}}))
    data = collect_report_data(tmp_path, path)
    assert data["warnings"] and data["screening"]["ranking"] == []

def test_provenance_paths_are_relative() -> None:
    assert all(not Path(path).is_absolute() for path in collect_report_data()["files"]["source_inputs"])
