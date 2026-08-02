"""Publication-oriented v2 text outputs."""
from __future__ import annotations
import json
from pathlib import Path
SECTIONS=["Methods","Dataset","Results","Limitations","Reproducibility checklist"]
def write_publication_bundle(output_dir,context):
 out=Path(output_dir); out.mkdir(parents=True,exist_ok=True); files=[]
 for section in SECTIONS:
  p=out/f"{section.lower().replace(' ','_')}.md"; p.write_text(f"# {section}\n\n{context.get(section,'Not evaluated.')}\n"); files.append(p)
 (out/"publication_manifest.json").write_text(json.dumps({"files":[p.name for p in files],"experimental_inactives_only":True},indent=2)+"\n"); return files
