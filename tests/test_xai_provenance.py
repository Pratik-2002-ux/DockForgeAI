import json
from pathlib import Path
from src.datasets.provenance import sha256_file, write_provenance

def test_checksum_is_deterministic(tmp_path: Path) -> None:
    path=tmp_path/"x"; path.write_text("abc"); assert sha256_file(path)==sha256_file(path) and len(sha256_file(path))==64
def test_provenance_manifest(tmp_path: Path) -> None:
    path=write_provenance(tmp_path/"p.json",[{"database":"test"}],"1"); assert json.loads(path.read_text())["source_release"]=="1"
