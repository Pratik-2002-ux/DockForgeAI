"""Cached validation artifact helpers for Streamlit."""
from __future__ import annotations
import io, zipfile
from pathlib import Path

def validation_bundle(root: Path) -> bytes:
    base=root/"data/processed/xai_validation"; stream=io.BytesIO()
    with zipfile.ZipFile(stream,"w",zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(base.rglob("*")):
            if path.is_file() and path.suffix.lower() in {".csv",".json",".md",".pdf",".html",".png",".svg"}: archive.write(path,path.relative_to(root))
    return stream.getvalue()

def external_template() -> bytes:
    return b"compound_id,canonical_smiles,true_label,source_database\nexample_1,CCO,,user_supplied\n"
