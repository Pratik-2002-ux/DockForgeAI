from pathlib import Path
import pytest
from pypdf import PdfReader
from src.reports.pdf_report import PDFReportError, generate_pdf_report
from src.reports.report_data import collect_report_data

def test_pdf_created_with_at_least_five_pages(tmp_path: Path) -> None:
    output = generate_pdf_report(collect_report_data(), tmp_path / "report.pdf")
    assert output.stat().st_size > 10_000 and len(PdfReader(output).pages) >= 5

def test_pdf_overwrite_protection(tmp_path: Path) -> None:
    output = tmp_path / "report.pdf"; output.write_bytes(b"existing")
    with pytest.raises(PDFReportError, match="exists"): generate_pdf_report(collect_report_data(), output)
