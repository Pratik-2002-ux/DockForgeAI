import io
import zipfile
from pathlib import Path
import pytest
from app.components.export_center import build_zip_bytes, format_file_size, safe_download_files

def test_file_size_formatting() -> None:
    assert format_file_size(0) == "0 B" and format_file_size(1536) == "1.5 KB"

def test_negative_file_size_rejected() -> None:
    with pytest.raises(ValueError): format_file_size(-1)

def test_download_filter_excludes_missing_and_external(tmp_path: Path) -> None:
    inside = tmp_path / "inside.txt"; inside.write_text("ok")
    assert safe_download_files([inside, tmp_path / "missing", Path("/etc/hosts")], tmp_path) == [inside]

def test_zip_archive_uses_relative_paths(tmp_path: Path) -> None:
    path = tmp_path / "dir/file.txt"; path.parent.mkdir(); path.write_text("ok")
    with zipfile.ZipFile(io.BytesIO(build_zip_bytes([path], tmp_path))) as archive:
        assert archive.namelist() == ["dir/file.txt"]
