from __future__ import annotations
import json
from pathlib import Path
def generate_model_cards(root:Path)->list[Path]:
 registry=json.loads((root/"data/processed/xai_models/model_registry.json").read_text())["models"];best=json.loads((root/"data/processed/xai_benchmarks/best_model.json").read_text());out=root/"data/processed/xai_model_cards";out.mkdir(parents=True,exist_ok=True);outputs=[]
 for entry in registry:
  payload={**entry,"split_method":"Bemis-Murcko scaffold group 70/15/15","intended_use":"Benchmark model-behaviour research","inappropriate_use":"Binding, therapeutic, clinical, or biological claims","class_labels":{"active":"ChEMBL potency-qualified","inactive_or_decoy":"DUD-E computational decoy, not experimentally proven inactive"},"limitations":["Source differences may create shortcuts","Good classification does not prove binding","Applicability-domain warnings must be considered"],"reproducibility_command":"conda run --no-capture-output -n dockforge python scripts/run_complete_xai_pipeline.py --offline"}
  path=out/f"{entry['model_id']}.json";path.write_text(json.dumps(payload,indent=2,default=str));(out/f"{entry['model_id']}.md").write_text(f"# {entry['model_id']}\n\nScaffold-safe benchmark model. DUD-E decoys are not experimentally proven inactive. Good classification does not prove binding or therapeutic activity.\n");outputs.append(path)
 return outputs
