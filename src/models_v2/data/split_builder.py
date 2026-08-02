"""Deterministic leakage-aware split assignments."""
from __future__ import annotations
import hashlib
import pandas as pd

def _group_split(values: pd.Series) -> pd.Series:
    unique=sorted(set(values.fillna("").astype(str))); scored=sorted(unique,key=lambda x:hashlib.sha256(("2026"+x).encode()).hexdigest()); n=len(scored); mapping={v:("train" if i<round(.7*n) else "validation" if i<round(.85*n) else "test") for i,v in enumerate(scored)}; return values.fillna("").astype(str).map(mapping)

def add_splits(frame:pd.DataFrame)->pd.DataFrame:
    out=frame.copy(); out["dataset_split_scaffold"]=_group_split(out.murcko_scaffold_smiles)
    source_groups=sorted(set(out.source_databases.astype(str)),key=lambda x:hashlib.sha256(("source-2026"+x).encode()).hexdigest())
    source_mapping={value:("test" if index==len(source_groups)-1 else "train") for index,value in enumerate(source_groups)}
    out["dataset_split_source"]=out.source_databases.astype(str).map(source_mapping); out["dataset_split_assay"]=_group_split(out.source_assay_ids)
    dates=pd.to_numeric(out.publication_date,errors="coerce"); order=dates.rank(method="first",pct=True); out["dataset_split_time"]=order.map(lambda x:"train" if pd.isna(x) or x<=.7 else "validation" if x<=.85 else "test"); return out

def leakage_count(frame:pd.DataFrame, group:str, split:str)->int:
    return int((frame.groupby(group)[split].nunique()>1).sum())
