"""Checksummed validation artifact manifest."""
from __future__ import annotations
from pathlib import Path
from .core import OUTPUT, ROOT, sha256, write_json

def build_manifest(root=OUTPUT):
    root=Path(root); files=[]
    relative_base = ROOT if root.is_relative_to(ROOT) else root
    for path in sorted(root.rglob("*")):
        if path.is_file() and path.name!="validation_manifest.json": files.append({"path":str(path.relative_to(relative_base)),"size_bytes":path.stat().st_size,"sha256":sha256(path)})
    manifest={"offline":True,"new_biological_data_downloaded":False,"locked_model_retrained":False,"artifacts":files}; write_json(root/"validation_manifest.json",manifest); return manifest
