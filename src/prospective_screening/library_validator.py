import pandas as pd
from rdkit import Chem
from .library_schema import REQUIRED_LIBRARY_COLUMNS
def validate_library(frame:pd.DataFrame):
 missing=REQUIRED_LIBRARY_COLUMNS-set(frame)
 if missing:raise ValueError(f"missing columns: {sorted(missing)}")
 if frame.compound_id.duplicated().any():raise ValueError("duplicate compound IDs")
 if frame.canonical_smiles.map(Chem.MolFromSmiles).isna().any():raise ValueError("invalid structure")
 return True
