"""Safe inference-only interface for future user-supplied data."""
from __future__ import annotations
import json
from pathlib import Path
import joblib, numpy as np, pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from scipy import sparse
from .core import BENCHMARKS, DATASET, FEATURES, ROOT, safe_metrics, write_json, sha256

def validate_external_dataset(input_path, id_column, smiles_column, output_dir, label_column=None, strict=False, allow_scaffold_overlap=False, overwrite=False):
    input_path=Path(input_path); output=Path(output_dir)
    if output.exists() and any(output.iterdir()) and not overwrite: raise FileExistsError(output)
    output.mkdir(parents=True,exist_ok=True); data=pd.read_parquet(input_path) if input_path.suffix.lower()==".parquet" else pd.read_csv(input_path)
    missing={id_column,smiles_column}-set(data.columns)
    if missing: raise ValueError(f"Missing required columns: {sorted(missing)}")
    locked=pd.read_csv(DATASET); exact=set(locked.canonical_isomeric_smiles); train=locked[locked.dataset_split.eq("train")]; train_scaff=set(train.murcko_scaffold_smiles.fillna("")); rows=[]; fps=[]
    from rdkit.Chem.Scaffolds import MurckoScaffold
    for _,r in data.iterrows():
        mol=Chem.MolFromSmiles(str(r[smiles_column])); valid=mol is not None; canonical=Chem.MolToSmiles(mol,isomericSmiles=True) if valid else None; scaffold=MurckoScaffold.MurckoScaffoldSmiles(mol=mol,includeChirality=True) if valid else None; rows.append({"compound_id":r[id_column],"canonical_smiles":canonical,"valid_structure":valid,"exact_training_overlap":canonical in exact if valid else False,"scaffold_training_overlap":scaffold in train_scaff if valid else False});
        if valid: fps.append(AllChem.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=True).GetFingerprint(mol))
    overlap=pd.DataFrame(rows); overlap.to_csv(output/"external_overlap_report.csv",index=False)
    if strict and (~overlap.valid_structure).any(): raise ValueError("Invalid structures in strict mode")
    usable=overlap.valid_structure & ~overlap.exact_training_overlap
    if not allow_scaffold_overlap: usable &= ~overlap.scaffold_training_overlap
    model_info=json.loads((BENCHMARKS/"best_model.json").read_text()); registry=json.loads((ROOT/"data/processed/xai_models/model_registry.json").read_text())["models"]; entry=next(x for x in registry if x["model_id"]==model_info["model_id"]); model=joblib.load(ROOT/entry["artifact_path"])
    valid_positions=np.flatnonzero(overlap.valid_structure); fp_by_pos=dict(zip(valid_positions,fps)); selected=np.flatnonzero(usable); matrix=np.zeros((len(selected),2048),dtype=np.uint8)
    for j,pos in enumerate(selected): DataStructs.ConvertToNumpyArray(fp_by_pos[pos],matrix[j])
    probabilities=model.predict_proba(sparse.csr_matrix(matrix))[:,1] if len(selected) else np.array([]); predictions=pd.DataFrame({"compound_id":data.iloc[selected][id_column].to_numpy(),"predicted_probability":probabilities,"predicted_class":np.where(probabilities>=.5,"active","inactive_or_decoy")})
    if label_column: predictions["true_label"]=data.iloc[selected][label_column].to_numpy()
    predictions.to_csv(output/"external_predictions.csv",index=False); ad=overlap.loc[selected,["compound_id","scaffold_training_overlap"]].copy(); ad["applicability_warning"]=np.where(ad.scaffold_training_overlap,"scaffold_overlap","not_assessed_beyond_overlap"); ad.to_csv(output/"external_applicability_domain.csv",index=False)
    metrics={"labels_supplied":bool(label_column),"model_retrained":False,"threshold_optimized":False,"usable_count":len(selected)}
    if label_column and len(selected): metrics.update(safe_metrics((predictions.true_label.astype(str)=="active").astype(int),probabilities))
    write_json(output/"external_validation_metrics.json",metrics); manifest={"input_name":input_path.name,"input_sha256":sha256(input_path),"locked_model_id":model_info["model_id"],"model_retrained":False,"locked_artifacts_overwritten":False,"independent_external_validation_complete":False,"outputs":[p.name for p in sorted(output.iterdir())]}; write_json(output/"external_validation_manifest.json",manifest); return manifest
