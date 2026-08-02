"""Canonical v2 dataset schema and validation."""
from __future__ import annotations
import pandas as pd
REQUIRED=["compound_id","canonical_isomeric_smiles","binary_label","source_database","source_assay_id","measurement_type","standard_value_nm","assay_date","murcko_scaffold_smiles","record_weight"]
EXPERIMENTAL_LABELS={"active","inactive"}
class V2SchemaError(ValueError): pass
def validate_dataset(frame:pd.DataFrame):
 missing=set(REQUIRED)-set(frame.columns)
 if missing: raise V2SchemaError(f"Missing v2 columns: {sorted(missing)}")
 if not set(frame.binary_label.dropna()).issubset(EXPERIMENTAL_LABELS): raise V2SchemaError("v2 accepts only experimentally measured active/inactive labels")
 if frame.canonical_isomeric_smiles.duplicated().any(): raise V2SchemaError("Duplicate standardized structures")
 if frame.standard_value_nm.le(0).any(): raise V2SchemaError("Activity values must be positive")
 return True
