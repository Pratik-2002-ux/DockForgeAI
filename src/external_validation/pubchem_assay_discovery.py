"""Discover and audit independent PubChem BioAssays."""
from __future__ import annotations
from datetime import datetime,timezone
import json
from pathlib import Path
import pandas as pd
from . import SEARCH_TERMS
from .pubchem_assay_auditor import audit_assay

def discover_assays(client,output_dir:Path,query="HIV-1 protease",exclude_depositor="ChEMBL",max_assays=250,minimum_usable=30,refresh=False):
    output=Path(output_dir); output.mkdir(parents=True,exist_ok=True); aids={162060}; searches=[]
    terms=[]
    for term in SEARCH_TERMS:
        terms.extend([term,f'"{term}" AND BindingDB',f'"{term}" AND ChEMBL'])
    if query not in terms: terms.insert(0,query)
    for term in terms:
        found,count,endpoint=client.search_assays(term,max_assays,refresh); aids.update(found); searches.append({"term":term,"reported_count":count,"retrieved_count":len(found),"endpoint":endpoint})
    summaries=client.assay_summaries(sorted(aids),refresh); audited=pd.DataFrame([audit_assay(r,minimum_usable) for r in summaries]).drop_duplicates("pubchem_aid").sort_values("pubchem_aid")
    accepted=audited[audited.assay_quality_status.ne("unsuitable")].copy(); rejected=audited[audited.assay_quality_status.eq("unsuitable")].copy(); audited.to_csv(output/"candidate_assays.csv",index=False); (output/"candidate_assays.json").write_text(audited.to_json(orient="records",indent=2)); accepted.to_csv(output/"accepted_assays.csv",index=False); rejected.to_csv(output/"rejected_assays.csv",index=False)
    manifest={"retrieval_timestamp":datetime.now(timezone.utc).isoformat(),"searches":searches,"candidate_count":len(audited),"accepted_count":len(accepted),"rejected_count":len(rejected),"explicit_exclusion":{"pubchem_aid":162060,"reason":"ChEMBL depositor; not independent"},"official_services_only":True}; (output/"assay_discovery_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n"); return manifest,audited
