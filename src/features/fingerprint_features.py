"""Morgan and MACCS sparse fingerprints."""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import sparse
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MACCSkeys

def fingerprint_matrix(dataset: pd.DataFrame, kind: str="morgan", n_bits: int=2048, radius: int=2, use_chirality: bool=True) -> sparse.csr_matrix:
    rows=[]; width=n_bits if kind=="morgan" else 167
    for smiles in dataset["canonical_isomeric_smiles"]:
        mol=Chem.MolFromSmiles(smiles)
        if mol is None: raise ValueError(f"Invalid standardized SMILES: {smiles}")
        fp=AllChem.GetMorganGenerator(radius=radius,fpSize=n_bits,includeChirality=use_chirality).GetFingerprint(mol) if kind=="morgan" else MACCSkeys.GenMACCSKeys(mol)
        array=np.zeros((width,),dtype=np.uint8); DataStructs.ConvertToNumpyArray(fp,array); rows.append(array)
    return sparse.csr_matrix(np.asarray(rows))
