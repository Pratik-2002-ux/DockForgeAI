"""Safe v2 artifact exports with explicit optional backends."""
from __future__ import annotations
import json,pickle
from pathlib import Path
def export_pickle(model,path): Path(path).write_bytes(pickle.dumps(model)); return Path(path)
def export_metadata(metadata,path): Path(path).write_text(json.dumps(metadata,indent=2)+"\n"); return Path(path)
def export_onnx(*args,**kwargs): raise RuntimeError("ONNX export requires the optional v2 ONNX backend")
def export_torchscript(*args,**kwargs): raise RuntimeError("TorchScript export requires the optional v2 PyTorch backend")
