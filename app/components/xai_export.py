from pathlib import Path
from app.components.export_center import build_zip_bytes
def xai_bundle(root:Path)->bytes:
 files=[]
 for directory in ("xai_predictions","xai_explanations","xai_benchmarks","xai_figures","xai_final_reports","xai_model_cards"):
  files.extend((root/"data/processed"/directory).glob("**/*"))
 return build_zip_bytes(files,root)
