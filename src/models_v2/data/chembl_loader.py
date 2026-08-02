"""Offline ChEMBL cache loader preserving record provenance."""
from __future__ import annotations
import json
from pathlib import Path
import pandas as pd
from .activity_harmonizer import classify_measurement, convert_to_nm, endpoint_fields
from .measurement_parser import parse_measurement
from .structure_standardizer import StructureStandardizationError, standardize_structure

BAD_VALIDITY=("outside typical range","invalid","ambiguous","conflict","transcription error")

def load_chembl(cache: str | Path, target="CHEMBL243") -> tuple[pd.DataFrame,pd.DataFrame]:
    records=[]
    for page in sorted(Path(cache).glob("activity_page_*.json")): records.extend(json.loads(page.read_text()).get("activities",[]))
    rows=[]; excluded=[]
    for idx,r in enumerate(records):
        reason=None; validity=str(r.get("data_validity_comment") or "").lower(); endpoint=str(r.get("standard_type") or "").upper()
        if r.get("target_chembl_id") != target: reason="wrong_target"
        elif any(x in validity for x in BAD_VALIDITY): reason="excluded_source_records"
        elif endpoint not in {"KI","KD","IC50","EC50"}: reason="incompatible_assays"
        relation,value=parse_measurement(r.get("standard_value"),r.get("standard_relation"))
        if reason is None and value is None: reason="missing_measurements"
        try: nm=convert_to_nm(value,r.get("standard_units")) if value is not None else None
        except ValueError: nm=None; reason="unsupported_units"
        try: structure=standardize_structure(str(r.get("canonical_smiles") or "")) if reason is None else {}
        except StructureStandardizationError: structure={}; reason="invalid_structures"
        if reason:
            excluded.append({"source_database":"ChEMBL","source_record_id":r.get("activity_id",idx),"exclusion_reason":reason,"detail":r.get("data_validity_comment")}); continue
        status=classify_measurement(nm,relation,r.get("activity_comment"),endpoint)
        rows.append({"source_database":"ChEMBL","source_record_id":str(r.get("activity_id") or idx),"source_compound_id":str(r.get("molecule_chembl_id") or ""),"source_assay_id":str(r.get("assay_chembl_id") or ""),"source_document_id":str(r.get("document_chembl_id") or ""),"source_publication_id":str(r.get("document_chembl_id") or ""),"source_depositor":"ChEMBL","target_name":str(r.get("target_pref_name") or ""),"target_identifier":target,"target_organism":str(r.get("target_organism") or ""),"target_variant":str(r.get("assay_variant_mutation") or ""),"assay_type":str(r.get("assay_type") or ""),"activity_type":endpoint,"original_value":r.get("standard_value"),"original_units":r.get("standard_units"),"activity_relation":relation,"standard_value_nm":nm,"pactivity":9-__import__('math').log10(nm),"categorical_outcome":str(r.get("activity_comment") or ""),"original_smiles":r.get("canonical_smiles"),**structure,"measurement_date":"","publication_date":r.get("document_year") or "","retrieval_date":"2026-08-02","provenance_reference":"data/raw/xai/chembl","source_quality_score":.9,"assay_quality_score":.9,"external_overlap_flag":False,"record_status":status,**endpoint_fields(endpoint,nm)})
    return pd.DataFrame(rows),pd.DataFrame(excluded)
