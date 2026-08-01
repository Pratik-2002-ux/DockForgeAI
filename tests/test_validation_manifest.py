from pathlib import Path
from src.validation.validation_manifest import build_manifest
def test_manifest_uses_relative_paths(tmp_path):
    (tmp_path/"x.csv").write_text("a\n1\n"); m=build_manifest(tmp_path); assert not Path(m["artifacts"][0]["path"]).is_absolute()
def test_manifest_confirms_offline(tmp_path): assert build_manifest(tmp_path)["offline"] is True
