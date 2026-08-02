"""Compound-level evidence aggregation without blind averaging."""
from __future__ import annotations
import pandas as pd

def resolve_compounds(measurements: pd.DataFrame) -> pd.DataFrame:
    rows=[]
    for key,g in measurements.groupby("inchikey_connectivity",sort=True):
        labels=set(g.record_status); values=g.pactivity.dropna(); active=int((g.record_status=="active").sum()); inactive=int((g.record_status=="inactive").sum()); ambiguous=int((g.record_status=="ambiguous").sum()); spread=float(values.max()-values.min()) if len(values) else float("nan")
        conflict=active>0 and inactive>0 or (pd.notna(spread) and spread>2)
        status="conflicting" if conflict else "active" if active else "inactive" if inactive else "ambiguous"
        representative=g.sort_values(["source_quality_score","assay_quality_score"],ascending=False).iloc[0]
        sets=lambda col:";".join(sorted(set(g[col].dropna().astype(str))-{""}))
        rows.append({"compound_id":f"V2_{key}","canonical_isomeric_smiles":representative.canonical_isomeric_smiles,"inchi":representative.inchi,"inchikey":representative.inchikey,"inchikey_connectivity":key,"murcko_scaffold_smiles":representative.murcko_scaffold_smiles,"binary_label":status if status in {"active","inactive"} else "","compound_status":status,"label_confidence":float(g.assay_quality_score.mean()),"measurement_count":len(g),"source_count":g.source_database.nunique(),"assay_count":g.source_assay_id.replace("",pd.NA).nunique(),"publication_count":g.source_publication_id.replace("",pd.NA).nunique(),"activity_types":sets("activity_type"),"minimum_pactivity":values.min() if len(values) else None,"maximum_pactivity":values.max() if len(values) else None,"median_pactivity":values.median() if len(values) else None,"pactivity_spread":spread,"primary_activity_type":representative.activity_type,"source_databases":sets("source_database"),"source_record_ids":sets("source_record_id"),"source_assay_ids":sets("source_assay_id"),"source_document_ids":sets("source_document_id"),"source_publication_ids":sets("source_publication_id"),"source_depositors":sets("source_depositor"),"assay_quality_score":float(g.assay_quality_score.mean()),"provenance_quality_score":float(g.source_quality_score.mean()),"active_record_count":active,"inactive_record_count":inactive,"ambiguous_record_count":ambiguous,"publication_date":pd.to_numeric(g.publication_date,errors="coerce").min()})
    return pd.DataFrame(rows)
