"""Exact Ki curation from the authenticated local BindingDB extract."""
from __future__ import annotations
import math
import pandas as pd
from src.models_v2.data.structure_standardizer import standardize_structure,StructureStandardizationError
from .ki_parser import parse_ki,pki_from_nm
def curate_bindingdb_ki(frame):
 exact=[]; censored=[]; invalid=[]
 for i,r in frame.iterrows():
  rel,nm=parse_ki(r.get("Ki (nM)"),"nM")
  if nm is None:
   if pd.notna(r.get("Ki (nM)")): invalid.append({"source_row":i,"reason":"malformed_or_nonpositive_ki"})
   continue
  try:s=standardize_structure(str(r.get("Ligand SMILES") or ""))
  except StructureStandardizationError as e: invalid.append({"source_row":i,"reason":"invalid_structure","detail":str(e)}); continue
  row={"source_row":i,"source_record_id":r.get("BindingDB Reactant_set_id"),"source_compound_id":r.get("BindingDB MonomerID"),"original_ki":r.get("Ki (nM)"),"ki_relation":rel,"ki_nm":nm,"pki":pki_from_nm(nm),"doi":r.get("Article DOI"),"institution":r.get("Institution"),"publication_year":r.get("Publication Year"),"target_name":r.get("Target Name"),**s}
  (exact if rel=="=" else censored).append(row)
 return pd.DataFrame(exact),pd.DataFrame(censored),pd.DataFrame(invalid)
def aggregate_exact(frame):
 rows=[]
 for key,g in frame.groupby("inchikey_connectivity",sort=True):
  vals=g.pki; rep=g.sort_values("source_record_id").iloc[0]; spread=float(vals.max()-vals.min())
  rows.append({"compound_id":f"KI_{key}","canonical_isomeric_smiles":rep.canonical_isomeric_smiles,"inchi":rep.inchi,"inchikey":rep.inchikey,"inchikey_connectivity":key,"murcko_scaffold_smiles":rep.murcko_scaffold_smiles,"minimum_ki_nm":g.ki_nm.min(),"maximum_ki_nm":g.ki_nm.max(),"median_ki_nm":g.ki_nm.median(),"median_pki":vals.median(),"pki_spread":spread,"measurement_count":len(g),"doi_count":g.doi.replace("",pd.NA).nunique(),"institution_count":g.institution.replace("",pd.NA).nunique(),"source_record_ids":";".join(sorted(g.source_record_id.astype(str))),"publication_year":pd.to_numeric(g.publication_year,errors="coerce").min(),"conflict_flag":spread>2})
 return pd.DataFrame(rows)
