import pandas as pd
from .library_schema import REQUIRED_LIBRARY_COLUMNS
def validate_library(frame:pd.DataFrame):
 missing=REQUIRED_LIBRARY_COLUMNS-set(frame)
 if missing:raise ValueError(f"missing columns: {sorted(missing)}")
 return True
