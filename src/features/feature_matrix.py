"""Build and persist descriptor, fingerprint, and fused matrices."""
from __future__ import annotations
import hashlib, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy import sparse
from src.features.descriptor_features import descriptor_frame, descriptor_quality
from src.features.fingerprint_features import fingerprint_matrix
from src.features.feature_schema import validate_feature_names

def build_feature_artifacts(dataset_path:Path,output_dir:Path,overwrite:bool=False)->dict:
    output_dir.mkdir(parents=True,exist_ok=True); manifest_path=output_dir/"feature_manifest.json"
    if manifest_path.exists() and not overwrite: raise FileExistsError(manifest_path)
    data=pd.read_csv(dataset_path); descriptors=descriptor_frame(data); quality=descriptor_quality(descriptors)
    constants=quality["constant_columns"]; descriptors=descriptors.drop(columns=constants); validate_feature_names(descriptors.columns.tolist())
    morgan=fingerprint_matrix(data); maccs=fingerprint_matrix(data,"maccs"); scaled=(descriptors-descriptors.mean())/descriptors.std(ddof=0).replace(0,1)
    fused=sparse.hstack([sparse.csr_matrix(scaled.to_numpy()),morgan,maccs],format="csr")
    descriptors.to_parquet(output_dir/"descriptor_features.parquet",index=False); sparse.save_npz(output_dir/"morgan_features.npz",morgan); sparse.save_npz(output_dir/"maccs_features.npz",maccs); sparse.save_npz(output_dir/"fused_features.npz",fused)
    index=data[["compound_id","standardized_compound_id","dataset_split","binary_label","scaffold_id"]]; index.to_csv(output_dir/"feature_index.csv",index=False)
    names={"descriptor":descriptors.columns.tolist(),"fingerprint":[f"Morgan_{i}" for i in range(morgan.shape[1])],"maccs":[f"MACCS_{i}" for i in range(maccs.shape[1])]}; schema={"feature_names":names,"dimensions":{"descriptor":list(descriptors.shape),"fingerprint":list(morgan.shape),"maccs":list(maccs.shape),"fused":list(fused.shape)}}
    (output_dir/"feature_schema.json").write_text(json.dumps(schema,indent=2)); (output_dir/"feature_quality_report.json").write_text(json.dumps(quality,indent=2))
    checksum=hashlib.sha256(dataset_path.read_bytes()).hexdigest(); manifest={"dataset_path":"data/processed/xai_dataset/"+dataset_path.name,"dataset_sha256":checksum,"row_count":len(data),"dimensions":schema["dimensions"],"constant_columns_removed":constants,"leakage_columns_included":[]}
    manifest_path.write_text(json.dumps(manifest,indent=2)); return manifest
