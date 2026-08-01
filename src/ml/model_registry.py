from __future__ import annotations
import json
from pathlib import Path
def write_registry(path:Path,entries:list[dict]): path.parent.mkdir(parents=True,exist_ok=True); path.write_text(json.dumps({"models":entries},indent=2)); return path
def load_registry(path:Path): return json.loads(path.read_text()).get("models",[])
