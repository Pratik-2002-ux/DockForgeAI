"""PubChem loader gate: locked external records are never training inputs."""
from pathlib import Path
import pandas as pd

def load_pubchem(path: str | Path | None, locked_root: str | Path) -> tuple[pd.DataFrame,pd.DataFrame]:
    if path is None: return pd.DataFrame(),pd.DataFrame()
    source=Path(path).resolve(); locked=Path(locked_root).resolve()
    if source == locked or locked in source.parents: raise ValueError("locked PubChem external data cannot be loaded for training")
    raise ValueError("No approved independent local PubChem training mapping is configured")
