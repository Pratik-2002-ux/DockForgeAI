"""Checksums and path-safe PubChem provenance manifests."""
from __future__ import annotations
import hashlib,json
from datetime import datetime,timezone
from pathlib import Path
def checksum(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write_manifest(path,files,extra=None):
    path=Path(path); root=Path.cwd(); entries=[]
    for f in files:
        f=Path(f); entries.append({"path":str(f.relative_to(root)) if f.is_relative_to(root) else f.name,"size_bytes":f.stat().st_size,"sha256":checksum(f)})
    value={"created_at":datetime.now(timezone.utc).isoformat(),"official_pubchem_services_only":True,"locked_model_retrained":False,"files":entries,**(extra or {})}; path.write_text(json.dumps(value,indent=2)+"\n"); return value
