import json
from pathlib import Path
import pytest
import yaml
from src.reports.final_report import FinalReportError, generate_final_reports

def _fixture_root(tmp_path: Path) -> Path:
    (tmp_path / "config").mkdir()
    keys = ("report_figures", "molecule_images", "pdf_reports", "html_reports", "final_reports")
    directories = {key: f"data/{key}" for key in keys}
    (tmp_path / "config/settings.yaml").write_text(yaml.safe_dump({"project": {"name": "DockForgeAI", "version": "test"}, "directories": directories}))
    return tmp_path

def test_manifest_created_with_temp_missing_fixtures(tmp_path: Path) -> None:
    root = _fixture_root(tmp_path)
    manifest = generate_final_reports(root, generate_pdf=False, generate_html=True)
    assert (root / manifest["manifest_path"]).is_file() and json.loads((root / manifest["manifest_path"]).read_text())["warnings"]

def test_strict_mode_rejects_missing_sections(tmp_path: Path) -> None:
    with pytest.raises(FinalReportError): generate_final_reports(_fixture_root(tmp_path), generate_pdf=False, strict=True)

def test_manifest_overwrite_protection(tmp_path: Path) -> None:
    root = _fixture_root(tmp_path); generate_final_reports(root, generate_pdf=False, generate_html=False)
    with pytest.raises(FinalReportError): generate_final_reports(root, generate_pdf=False, generate_html=False)
