#!/usr/bin/env python3
from __future__ import annotations
import argparse,hashlib,json,sys
from pathlib import Path
import pandas as pd
from rdkit import Chem,RDLogger
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from src.models_v2.regression.ki_curator import curate_bindingdb_ki,aggregate_exact
from src.models_v2.regression.external_overlap import locked_identities,mark_external_overlap
from src.models_v2.regression.scaffold_splits import assign_scaffold_split
def main():
 p=argparse.ArgumentParser();p.add_argument("--input",default="data/model_v2/raw/bindingdb_hiv1_protease_filtered.csv");p.add_argument("--locked",default="data/processed/pubchem_external_validation/dataset/pubchem_external_validation.csv");p.add_argument("--output-dir",default="data/model_v2/regression/dataset");p.add_argument("--overwrite",action="store_true");a=p.parse_args();RDLogger.DisableLog("rdApp.*");out=ROOT/a.output_dir;out.mkdir(parents=True,exist_ok=True)
 raw=pd.read_csv(ROOT/a.input,dtype=str,low_memory=False);exact,censored,invalid=curate_bindingdb_ki(raw);comp=aggregate_exact(exact);comp=mark_external_overlap(comp,locked_identities(ROOT/a.locked));conf=comp[comp.conflict_flag];over=comp[comp.external_overlap_flag];final=assign_scaffold_split(comp[~comp.conflict_flag&~comp.external_overlap_flag].copy()).sort_values("compound_id").reset_index(drop=True)
 files={"dataset":out/"ki_regression_dataset.csv","parquet":out/"ki_regression_dataset.parquet","sdf":out/"ki_regression_dataset.sdf","measurements":out/"ki_measurement_table.csv","censored":out/"ki_censored_measurements.csv","conflicts":out/"ki_conflicting_measurements.csv","invalid":out/"ki_invalid_structures.csv","overlap":out/"ki_locked_overlap.csv"};final.to_csv(files["dataset"],index=False);final.to_parquet(files["parquet"],index=False);exact.to_csv(files["measurements"],index=False);censored.to_csv(files["censored"],index=False);conf.to_csv(files["conflicts"],index=False);invalid.to_csv(files["invalid"],index=False);over.to_csv(files["overlap"],index=False)
 w=Chem.SDWriter(str(files["sdf"]));
 for _,r in final.iterrows():
  m=Chem.MolFromSmiles(r.canonical_isomeric_smiles)
  for k,v in r.items():m.SetProp(str(k),"" if pd.isna(v) else str(v))
  w.write(m)
 w.close();stats={"raw_rows":len(raw),"exact_measurements":len(exact),"censored_measurements":len(censored),"invalid_measurements":len(invalid),"unique_before_exclusion":len(comp),"conflicting_compounds":len(conf),"locked_overlap_compounds":len(over),"final_compounds":len(final),"unique_scaffolds":final.murcko_scaffold_smiles.nunique(),"split_counts":final.dataset_split.value_counts().to_dict(),"scaffold_leakage":int((final.groupby('murcko_scaffold_smiles').dataset_split.nunique()>1).sum())}
 manifest={"endpoint":"Ki","ordinary_regression_relations":["="],"censored_excluded":True,"locked_external_training_use":False,"stats":stats};(out/"ki_dataset_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n");(out/"ki_dataset_schema.json").write_text(json.dumps({"columns":{c:str(t) for c,t in final.dtypes.items()}},indent=2)+"\n");checks={str(x.relative_to(ROOT)):hashlib.sha256(x.read_bytes()).hexdigest() for x in files.values()};(out/"ki_dataset_checksums.json").write_text(json.dumps(checks,indent=2)+"\n");print(json.dumps(stats,indent=2));return 0
if __name__=="__main__":raise SystemExit(main())
