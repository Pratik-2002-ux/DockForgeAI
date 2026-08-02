"""PubChem external-validation dashboard exports."""
from __future__ import annotations
import io,zipfile
from pathlib import Path
DOWNLOAD_KEYS={"pdf":"pubchem_external_pdf","html":"pubchem_external_html","dataset":"pubchem_external_dataset","bundle":"pubchem_external_bundle"}
def build_bundle(root:Path):
    stream=io.BytesIO(); base=root/"data/processed/pubchem_external_validation"
    with zipfile.ZipFile(stream,"w",zipfile.ZIP_DEFLATED) as z:
        for p in sorted(base.rglob("*")):
            if p.is_file() and p.suffix.lower() in {".csv",".json",".pdf",".html"}: z.write(p,p.relative_to(root))
    return stream.getvalue()
