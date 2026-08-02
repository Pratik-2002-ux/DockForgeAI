#!/usr/bin/env python3
"""Build the offline v2 multi-source experimental dataset; never trains."""
from __future__ import annotations
import argparse,hashlib,json,sys
from collections import Counter
from pathlib import Path
import pandas as pd
from rdkit import Chem, RDLogger

ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))
from src.models_v2.data.bindingdb_loader import load_bindingdb
from src.models_v2.data.chembl_loader import load_chembl
from src.models_v2.data.dataset_merger import merge_measurements
from src.models_v2.data.conflict_resolver import resolve_compounds
from src.models_v2.data.external_overlap import locked_identities,mark_external_overlap
from src.models_v2.data.split_builder import add_splits
from src.models_v2.data.dataset_report import readiness

EXCLUSIONS=["invalid_structures","missing_measurements","unsupported_units","ambiguous_measurements","conflicting_measurements","exact_external_overlap","connectivity_external_overlap","duplicate_structures","missing_provenance","incompatible_assays","excluded_source_records"]
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def write_sdf(path,frame):
    writer=Chem.SDWriter(str(path))
    for _,r in frame.iterrows():
        mol=Chem.MolFromSmiles(r.canonical_isomeric_smiles)
        if mol:
            for key,value in r.items(): mol.SetProp(str(key),"" if pd.isna(value) else str(value))
            writer.write(mol)
    writer.close()
def main():
    RDLogger.DisableLog("rdApp.*")
    p=argparse.ArgumentParser(); p.add_argument("--bindingdb",required=True); p.add_argument("--chembl-cache",required=True); p.add_argument("--locked-external",required=True); p.add_argument("--output-dir",required=True); p.add_argument("--overwrite",action="store_true"); a=p.parse_args(); out=ROOT/a.output_dir; out.mkdir(parents=True,exist_ok=True); (out/"exclusions").mkdir(exist_ok=True)
    b,bx=load_bindingdb(ROOT/a.bindingdb); c,cx=load_chembl(ROOT/a.chembl_cache); measurements=merge_measurements(b,c); measurements=mark_external_overlap(measurements,locked_identities(ROOT/a.locked_external)); compounds=resolve_compounds(measurements); compounds=mark_external_overlap(compounds,locked_identities(ROOT/a.locked_external))
    exact=compounds[compounds.external_exact_smiles_overlap|compounds.external_inchikey_overlap]; connectivity=compounds[~(compounds.external_exact_smiles_overlap|compounds.external_inchikey_overlap)&compounds.external_connectivity_overlap]; ambiguous=compounds[compounds.compound_status=="ambiguous"]; conflicts=compounds[compounds.compound_status=="conflicting"]
    final=compounds[(compounds.compound_status.isin(["active","inactive"]))&~compounds.external_overlap_flag].copy(); final["locked_external_overlap"]=False; final=add_splits(final).sort_values("compound_id").reset_index(drop=True)
    paths={"csv":out/"v2_multisource_experimental_dataset.csv","parquet":out/"v2_multisource_experimental_dataset.parquet","sdf":out/"v2_multisource_experimental_dataset.sdf","measurements":out/"v2_multisource_measurement_table.csv","summary":out/"v2_multisource_compound_summary.csv"}
    final.to_csv(paths["csv"],index=False); final.to_parquet(paths["parquet"],index=False); write_sdf(paths["sdf"],final); measurements.to_csv(paths["measurements"],index=False); compounds.to_csv(paths["summary"],index=False)
    exclusions=pd.concat([bx,cx],ignore_index=True); mapping={name:exclusions[exclusions.exclusion_reason==name] if not exclusions.empty else pd.DataFrame() for name in EXCLUSIONS}; mapping["ambiguous_measurements"]=ambiguous; mapping["conflicting_measurements"]=conflicts; mapping["exact_external_overlap"]=exact; mapping["connectivity_external_overlap"]=connectivity; mapping["duplicate_structures"]=measurements[measurements.duplicated(["source_database","source_record_id","activity_type","standard_value_nm"],keep=False)]
    for name,frame in mapping.items(): frame.to_csv(out/"exclusions"/f"{name}.csv",index=False)
    by_source={s:{"measurement_rows":len(g),"active":int((g.record_status=="active").sum()),"inactive":int((g.record_status=="inactive").sum()),"ambiguous":int((g.record_status=="ambiguous").sum()),"unique_structures":g.inchikey_connectivity.nunique(),"assays":g.source_assay_id.nunique(),"publications":g.source_publication_id.replace("",pd.NA).nunique()} for s,g in measurements.groupby("source_database")}
    summary={"bindingdb_rows_used":len(b),"chembl_rows_used":len(c),"pubchem_rows_used":0,"source_counts":by_source,"final_active_count":int((final.binary_label=="active").sum()),"final_inactive_count":int((final.binary_label=="inactive").sum()),"ambiguous_count":len(ambiguous),"conflicting_count":len(conflicts),"unique_compound_count":len(final),"unique_scaffold_count":final.murcko_scaffold_smiles.nunique(),"source_count":measurements.source_database.nunique(),"assay_count":measurements.source_assay_id.nunique(),"publication_count":measurements.source_publication_id.replace("",pd.NA).nunique(),"exact_external_overlap_count":len(exact),"connectivity_external_overlap_count":len(connectivity),"scaffold_external_overlap_count":int(final.locked_external_scaffold_overlap.sum()),"locked_external_overlap_count":int(final.external_overlap_flag.sum()),"duplicate_measurement_count":len(mapping["duplicate_structures"]),"split_counts":{col:dict(Counter(final[col])) for col in ["dataset_split_scaffold","dataset_split_source","dataset_split_assay","dataset_split_time"]}}
    source_label=pd.crosstab(final.source_databases,final.binary_label); source_label.to_csv(out/"source_label_table.csv"); shortcut=any((source_label.get(label,pd.Series(dtype=float))/max(1,source_label.get(label,pd.Series(dtype=float)).sum())>.8).any() for label in ["active","inactive"]); classical,gnn=readiness(summary); summary.update({"source_shortcut_risk":bool(shortcut),"classical_training_ready":classical and not shortcut,"gnn_training_ready":gnn and not shortcut}); summary["blockers"]=[]
    if summary["final_inactive_count"]<200: summary["blockers"].append(f"Measured inactive deficit: {200-summary['final_inactive_count']}")
    if shortcut: summary["blockers"].append("One source contributes more than 80% of a class.")
    summary["training_allowed"]=bool(summary["classical_training_ready"])
    summary["status"]="dataset_ready" if summary["training_allowed"] else "readiness_only"
    summary["approved_training_sources"]=["BindingDB","ChEMBL"]
    summary["locked_external"]={"path":a.locked_external,"exists":True,"training_use_forbidden":True}
    summary["warnings"]=["The locked PubChem/BindingDB benchmark was exclusion-only and never merged as training data."]; summary["required_next_action"]="Do not train; obtain approved independent measured inactive records." if summary["blockers"] else "Dataset is ready for a separately authorized modelling milestone."
    readiness_path=ROOT/"data/model_v2/readiness/readiness_report.json"; readiness_path.parent.mkdir(parents=True,exist_ok=True); readiness_path.write_text(json.dumps(summary,indent=2,default=int)+"\n")
    manifest={"schema_version":1,"offline":True,"training_performed":False,"locked_external_training_use":False,"summary":summary,"inputs":{"bindingdb":a.bindingdb,"chembl_cache":a.chembl_cache,"locked_external":a.locked_external}}
    (out/"v2_multisource_dataset_manifest.json").write_text(json.dumps(manifest,indent=2,default=int)+"\n"); (out/"v2_multisource_dataset_schema.json").write_text(json.dumps({"columns":{c:str(t) for c,t in final.dtypes.items()}},indent=2)+"\n"); checks={str(p.relative_to(ROOT)):sha(p) for p in paths.values()}; (out/"v2_multisource_dataset_checksums.json").write_text(json.dumps(checks,indent=2)+"\n")
    audit=ROOT/"data/model_v2/audit"; audit.mkdir(parents=True,exist_ok=True)
    for name,data in [("bindingdb_source_audit",by_source.get("BindingDB",{})),("chembl_source_audit",by_source.get("ChEMBL",{})),("pubchem_source_audit",{"rows_used":0,"reason":"No approved independent local PubChem training records."}),("source_independence_audit",{"locked_external_training_use":False,"exact_removed":len(exact),"connectivity_removed":len(connectivity)}),("source_column_mapping",{"BindingDB":"inspected 633-column export","ChEMBL":"cached activities JSON common schema"})]: (audit/f"{name}.json").write_text(json.dumps(data,indent=2,default=int)+"\n")
    print(json.dumps(summary,indent=2,default=int)); return 0
if __name__=="__main__": raise SystemExit(main())
