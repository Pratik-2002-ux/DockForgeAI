"""Normalize locally supplied ChEMBL, PubChem, and BindingDB exports."""
from __future__ import annotations
from pathlib import Path
import pandas as pd
SUPPORTED_SOURCES={"ChEMBL","PubChem BioAssay","BindingDB"}
def load_source(path,source,column_map):
 if source not in SUPPORTED_SOURCES: raise ValueError(f"Unsupported v2 source: {source}")
 path=Path(path); frame=pd.read_parquet(path) if path.suffix.lower()==".parquet" else pd.read_csv(path); missing=set(column_map)-set(frame.columns)
 if missing: raise ValueError(f"Missing source columns: {sorted(missing)}")
 result=frame.rename(columns=column_map).copy(); result["source_database"]=source; result["source_file"]=path.name; return result
def integrate_sources(items): return pd.concat([load_source(path,source,mapping) for path,source,mapping in items],ignore_index=True)
