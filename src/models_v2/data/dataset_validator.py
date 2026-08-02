"""Dataset consistency and leakage validation."""
from __future__ import annotations
from pathlib import Path
import pandas as pd
from rdkit import Chem
from .split_builder import leakage_count

def validate_dataset(path:str|Path)->dict[str,object]:
    p=Path(path); frame=pd.read_csv(p); parquet=pd.read_parquet(p.with_suffix(".parquet")); sdf=list(Chem.SDMolSupplier(str(p.with_suffix(".sdf")),removeHs=False)); result={"row_count":len(frame),"parquet_row_count":len(parquet),"sdf_row_count":sum(m is not None for m in sdf),"duplicate_structures":int(frame.canonical_isomeric_smiles.duplicated().sum()),"locked_external_overlap_count":int(frame.locked_external_overlap.sum()),"scaffold_leakage_count":leakage_count(frame,"murcko_scaffold_smiles","dataset_split_scaffold"),"source_leakage_count":leakage_count(frame,"source_databases","dataset_split_source"),"assay_leakage_count":leakage_count(frame,"source_assay_ids","dataset_split_assay")}
    if len({result["row_count"],result["parquet_row_count"],result["sdf_row_count"]})>1 or any(result[k] for k in ("duplicate_structures","locked_external_overlap_count","scaffold_leakage_count","source_leakage_count","assay_leakage_count")): raise ValueError(result)
    return result
