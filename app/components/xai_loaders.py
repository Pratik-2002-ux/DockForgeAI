from __future__ import annotations
import json
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]
def load_xai_json(relative:str):
 path=ROOT/relative
 try:return json.loads(path.read_text())
 except (OSError,json.JSONDecodeError):return {}
def load_xai_csv(relative:str):
 try:return pd.read_csv(ROOT/relative)
 except (OSError,pd.errors.ParserError):return pd.DataFrame()
