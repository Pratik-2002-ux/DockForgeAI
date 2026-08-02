"""BindingDB mapping using the inspected local export schema."""
from __future__ import annotations
from pathlib import Path
import pandas as pd
from .activity_harmonizer import classify_measurement, endpoint_fields
from .measurement_parser import parse_measurement
from .structure_standardizer import StructureStandardizationError, standardize_structure

ENDPOINTS = {"Ki (nM)":"KI", "Kd (nM)":"KD", "IC50 (nM)":"IC50", "EC50 (nM)":"EC50"}

def load_bindingdb(path: str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw=pd.read_csv(path,dtype=str,low_memory=False); rows=[]; excluded=[]
    for idx,r in raw.iterrows():
        try: structure=standardize_structure(str(r.get("Ligand SMILES") or ""))
        except StructureStandardizationError as exc:
            excluded.append({"source_database":"BindingDB","source_row":idx,"exclusion_reason":"invalid_structures","detail":str(exc)}); continue
        found=False
        for column,endpoint in ENDPOINTS.items():
            relation,value=parse_measurement(r.get(column))
            if value is None: continue
            found=True; status=classify_measurement(value,relation,endpoint=endpoint)
            rows.append({"source_database":"BindingDB","source_record_id":str(r.get("BindingDB Reactant_set_id") or idx),"source_compound_id":str(r.get("BindingDB MonomerID") or ""),"source_assay_id":str(r.get("BindingDB Reactant_set_id") or ""),"source_document_id":str(r.get("BindingDB Entry DOI") or r.get("Article DOI") or ""),"source_publication_id":str(r.get("Article DOI") or ""),"source_depositor":str(r.get("Institution") or "BindingDB"),"target_name":str(r.get("Target Name") or ""),"target_identifier":"HIV-1 protease","target_organism":str(r.get("Target Source Organism According to Curator or DataSource") or ""),"target_variant":str(r.get("_possible_mutant_or_variant") or ""),"assay_type":"binding","activity_type":endpoint,"original_value":r.get(column),"original_units":"nM","activity_relation":relation,"standard_value_nm":value,"pactivity":9-__import__('math').log10(value),"categorical_outcome":"","original_smiles":r.get("Ligand SMILES"),**structure,"measurement_date":"","publication_date":"","retrieval_date":"2026-08-02","provenance_reference":"data/model_v2/raw/bindingdb_hiv1_protease_filtered.csv","source_quality_score":.9,"assay_quality_score":.85,"external_overlap_flag":False,"record_status":status,**endpoint_fields(endpoint,value)})
        if not found: excluded.append({"source_database":"BindingDB","source_row":idx,"exclusion_reason":"missing_measurements"})
    return pd.DataFrame(rows),pd.DataFrame(excluded)
