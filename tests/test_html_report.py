from pathlib import Path
import pytest
from src.reports.html_report import HTMLReportError, generate_html_report
from src.reports.report_data import collect_report_data

def test_html_has_required_sections_and_no_private_path(tmp_path: Path) -> None:
    text = generate_html_report(collect_report_data(), tmp_path / "report.html").read_text()
    assert "Executive summary" in text and "Scientific limitations" in text and "Reproducibility" in text
    assert "/Users/" not in text and "data/processed/" in text

def test_html_embeds_image(tmp_path: Path) -> None:
    image = tmp_path / "test.png"; image.write_bytes(b"png-data")
    assert "data:image/png;base64" in generate_html_report(collect_report_data(), tmp_path / "report.html", [image]).read_text()

def test_html_overwrite_protection(tmp_path: Path) -> None:
    output = tmp_path / "report.html"; output.write_text("existing")
    with pytest.raises(HTMLReportError): generate_html_report(collect_report_data(), output)
