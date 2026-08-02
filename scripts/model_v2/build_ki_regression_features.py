#!/usr/bin/env python3
import json,sys
from pathlib import Path
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from src.models_v2.regression.regression_features import matrices,DESCRIPTORS
def main():
 d=pd.read_csv(ROOT/"data/model_v2/regression/dataset/ki_regression_dataset.csv");out=ROOT/"data/model_v2/regression/features";out.mkdir(parents=True,exist_ok=True);m=matrices(d.canonical_isomeric_smiles);pd.DataFrame(m["descriptors"],columns=[x[0] for x in DESCRIPTORS]).to_parquet(out/"descriptors.parquet",index=False)
 for k,v in m.items():
  if k!="descriptors":np.savez_compressed(out/(k+".npz"),X=v)
 d[["compound_id","dataset_split"]].to_csv(out/"feature_index.csv",index=False);(out/"feature_schema.json").write_text(json.dumps({k:list(v.shape) for k,v in m.items()},indent=2)+"\n");(out/"feature_manifest.json").write_text(json.dumps({"target_leakage":False,"provenance_features":False,"representations":list(m)},indent=2)+"\n");print({k:v.shape for k,v in m.items()});return 0
if __name__=="__main__":raise SystemExit(main())
