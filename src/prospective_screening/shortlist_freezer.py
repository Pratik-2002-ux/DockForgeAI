from __future__ import annotations
import hashlib,json
from pathlib import Path
def freeze_shortlist(frame,path,metadata):
 path=Path(path);path.mkdir(parents=True,exist_ok=True);csv=path/"shortlisted_candidates.csv"
 if csv.exists():raise FileExistsError("shortlist is immutable")
 frame.sort_values("compound_id").to_csv(csv,index=False);checksum=hashlib.sha256(csv.read_bytes()).hexdigest();manifest={**metadata,"shortlist_checksum":checksum,"frozen":True,"experimental_results_seen":False};(path/"shortlist_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n");return manifest
