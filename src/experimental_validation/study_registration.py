from __future__ import annotations
import json
from pathlib import Path
REQUIRED={"study_id","shortlist_checksum","laboratory_id","assay_sop_id","responsible_investigator","hit_confirmation_rule","blinding_status","analysis_version"}
def register(data,path):
 missing=REQUIRED-set(data)
 if missing:raise ValueError(f"missing registration fields: {sorted(missing)}")
 path=Path(path)
 if path.exists():raise FileExistsError("registration is immutable")
 path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(data,indent=2)+"\n");return data
