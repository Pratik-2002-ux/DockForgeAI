from pathlib import Path
from app.components.xai_export import xai_bundle
def test_empty_xai_bundle_valid(tmp_path:Path): assert xai_bundle(tmp_path).startswith(b"PK")
