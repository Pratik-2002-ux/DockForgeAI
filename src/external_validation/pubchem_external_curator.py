"""Curate accepted PubChem assays into an overlap-audited external dataset."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np,pandas as pd
from rdkit import Chem,DataStructs
from rdkit.Chem import AllChem
from .pubchem_structure_resolver import resolve_result_structures
from .pubchem_provenance import write_manifest,checksum

OUTCOME_MAP={"active":"active","inactive":"inactive_or_decoy"}
def map_outcome(value): return OUTCOME_MAP.get(str(value).strip().lower())
def build_external_dataset(client,accepted_path:Path,output_dir:Path,training_path:Path,minimum_usable=30,refresh=False):
    output=Path(output_dir); output.mkdir(parents=True,exist_ok=True); exclusions=output/"exclusions"; exclusions.mkdir(exist_ok=True); accepted=pd.read_csv(accepted_path); training=pd.read_csv(training_path); exact_smiles=set(training.canonical_isomeric_smiles.dropna()); exact_keys=set(training.inchikey.dropna()); connectivity=set(training.inchikey_connectivity.dropna()); scaffolds=set(training.murcko_scaffold_smiles.dropna()); generator=AllChem.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=True); train_fps=[generator.GetFingerprint(Chem.MolFromSmiles(x)) for x in training.canonical_isomeric_smiles]; curated=[]; excluded={k:[] for k in ("inconclusive_outcomes","unspecified_outcomes","missing_structures","invalid_structures","duplicate_structures","contradictory_measurements","exact_training_overlap","ambiguous_assay_annotations","unsupported_outcome_types")}
    for _,assay in accepted.iterrows():
        aid=int(assay.pubchem_aid); raw=client.assay_results(aid,refresh); frame=resolve_result_structures(client,raw,refresh); seen=set()
        for _,row in frame.iterrows():
            outcome=str(row.get("PUBCHEM_ACTIVITY_OUTCOME","")).strip(); label=map_outcome(outcome); base={"pubchem_aid":aid,"pubchem_sid":row.get("PUBCHEM_SID"),"pubchem_cid":row.get("PUBCHEM_CID"),"pubchem_outcome":outcome}
            if outcome.lower()=="inconclusive": excluded["inconclusive_outcomes"].append(base); continue
            if outcome.lower() in {"unspecified",""}: excluded["unspecified_outcomes"].append(base); continue
            if label is None: excluded["unsupported_outcome_types"].append(base); continue
            if not row.get("structure_valid"): excluded["missing_structures" if not row.get("resolved_smiles") else "invalid_structures"].append(base); continue
            smiles=row.canonical_isomeric_smiles; key=row.inchikey
            if smiles in seen: excluded["duplicate_structures"].append(base); continue
            seen.add(smiles); exact=smiles in exact_smiles or key in exact_keys
            if exact: excluded["exact_training_overlap"].append({**base,"canonical_smiles":smiles}); continue
            mol=Chem.MolFromSmiles(smiles); fp=generator.GetFingerprint(mol); similarities=sorted(DataStructs.BulkTanimotoSimilarity(fp,train_fps),reverse=True); scaffold=row.murcko_scaffold_smiles; conn=row.inchikey_connectivity in connectivity
            curated.append({"pubchem_aid":aid,"pubchem_sid":row.get("PUBCHEM_SID"),"pubchem_cid":row.get("PUBCHEM_CID"),"compound_id":f"PUBCHEM_AID{aid}_SID{row.get('PUBCHEM_SID')}","assay_name":assay.assay_name,"depositor":assay.depositor,"target_name":assay.target_name,"assay_type":assay.assay_type,"pubchem_outcome":outcome,"pubchem_activity_score":row.get("PUBCHEM_ACTIVITY_SCORE"),"true_label":label,"original_measurement":row.get("PubChem Standard Value"),"original_units":"uM" if pd.notna(row.get("PubChem Standard Value")) else None,"original_qualifier":row.get("Standard Relation"),"canonical_smiles":smiles,"canonical_isomeric_smiles":smiles,"inchi":row.inchi,"inchikey":key,"inchikey_connectivity":row.inchikey_connectivity,"murcko_scaffold_smiles":scaffold,**{c:row.get(c) for c in ("molecular_weight","exact_molecular_weight","logp","tpsa","hbd","hba","rotatable_bonds","ring_count","aromatic_ring_count","formal_charge","heavy_atom_count","fraction_csp3","qed","pains_alert_count","brenk_alert_count")},"exact_training_overlap":False,"connectivity_training_overlap":conn,"scaffold_training_overlap":scaffold in scaffolds if scaffold else False,"maximum_training_tanimoto":similarities[0],"mean_top5_training_tanimoto":float(np.mean(similarities[:5])),"applicability_domain_status":"in_domain" if similarities[0]>=.4 else "out_of_domain","independence_class":"possible_provenance_overlap" if conn else ("structure-independent but scaffold-overlapping" if scaffold in scaffolds else "novel-scaffold"),"source_publication":None,"source_endpoint":f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/{aid}/CSV","retrieval_date":pd.Timestamp.utcnow().date().isoformat()})
    data=pd.DataFrame(curated); counts=data.groupby("pubchem_aid").size() if not data.empty else pd.Series(dtype=int); class_counts=data.groupby("pubchem_aid").true_label.nunique() if not data.empty else pd.Series(dtype=int); insufficient=set(counts[counts<minimum_usable].index)|set(class_counts[class_counts<2].index)
    audit_rows=[]
    for _,assay in accepted.iterrows():
        aid=int(assay.pubchem_aid); subset=data[data.pubchem_aid.eq(aid)]; exact_count=sum(int(r.get("pubchem_aid",-1))==aid for r in excluded["exact_training_overlap"]); reasons=[]
        if len(subset)<minimum_usable: reasons.append(f"fewer_than_{minimum_usable}_usable_compounds_after_overlap_removal")
        if subset.true_label.nunique()<2: reasons.append("no_inactive_or_active_outcomes_after_curation")
        row=assay.to_dict(); row.update(structure_mapping_available=bool(len(subset)),training_exact_overlap_count=exact_count,training_scaffold_overlap_count=int(subset.scaffold_training_overlap.sum()) if len(subset) else 0,usable_compound_count=len(subset),assay_quality_status="unsuitable" if reasons else "acceptable_with_caution",rejection_reasons="; ".join(reasons)); audit_rows.append(row)
    final_audit=pd.DataFrame(audit_rows); final_accepted=final_audit[final_audit.assay_quality_status.ne("unsuitable")]; newly_rejected=final_audit[final_audit.assay_quality_status.eq("unsuitable")]; final_accepted.to_csv(accepted_path,index=False); rejected_path=Path(accepted_path).with_name("rejected_assays.csv"); previous=pd.read_csv(rejected_path) if rejected_path.exists() else pd.DataFrame(); pd.concat([previous,newly_rejected],ignore_index=True).drop_duplicates("pubchem_aid",keep="last").to_csv(rejected_path,index=False)
    if insufficient: data=data[~data.pubchem_aid.isin(insufficient)]
    csv=output/"pubchem_external_validation.csv"; parquet=output/"pubchem_external_validation.parquet"; sdf=output/"pubchem_external_validation.sdf"; data.to_csv(csv,index=False); data.to_parquet(parquet,index=False)
    writer=Chem.SDWriter(str(sdf))
    for _,row in data.iterrows():
        mol=Chem.MolFromSmiles(row.canonical_isomeric_smiles); mol.SetProp("compound_id",row.compound_id); mol.SetProp("true_label",row.true_label); writer.write(mol)
    writer.close()
    for name,rows in excluded.items(): pd.DataFrame(rows).to_csv(exclusions/f"{name}.csv",index=False)
    schema={"columns":{c:str(t) for c,t in data.dtypes.items()},"required_columns":list(data.columns)}; (output/"pubchem_external_dataset_schema.json").write_text(json.dumps(schema,indent=2)+"\n"); checks={p.name:checksum(p) for p in (csv,parquet,sdf)}; (output/"pubchem_external_dataset_checksums.json").write_text(json.dumps(checks,indent=2)+"\n"); manifest=write_manifest(output/"pubchem_external_dataset_manifest.json",[csv,parquet,sdf],{"row_count":len(data),"accepted_aids":sorted(data.pubchem_aid.unique().tolist()) if not data.empty else [],"insufficient_after_curation_aids":sorted(map(int,insufficient)),"model_retrained":False}); return manifest
