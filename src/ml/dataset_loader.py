from __future__ import annotations
from pathlib import Path
import pandas as pd
from scipy import sparse
def load_feature_set(feature_dir:Path,name:str):
 index=pd.read_csv(feature_dir/"feature_index.csv"); y=(index.binary_label=="active").astype(int).to_numpy()
 if name=="descriptor": x=pd.read_parquet(feature_dir/"descriptor_features.parquet").to_numpy()
 elif name=="fingerprint": x=sparse.load_npz(feature_dir/"morgan_features.npz")
 elif name=="fused": x=sparse.load_npz(feature_dir/"fused_features.npz")
 else: raise ValueError(name)
 if x.shape[0]!=len(index): raise ValueError("Feature/index row mismatch")
 return x,y,index
def split_indices(index:pd.DataFrame)->dict[str,list[int]]: return {name:index.index[index.dataset_split==name].tolist() for name in ("train","validation","test")}
