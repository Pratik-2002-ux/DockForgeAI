"""Resolve and standardize PubChem assay structures."""
from __future__ import annotations
import pandas as pd
from src.datasets.structure_standardizer import standardize_structure, StructureStandardizationError

def resolve_result_structures(client,results:pd.DataFrame,refresh=False):
    frame=results.copy(); cid_col="PUBCHEM_CID"; smiles_col="PUBCHEM_EXT_DATASOURCE_SMILES"; cids=pd.to_numeric(frame.get(cid_col),errors="coerce").dropna().astype(int).unique().tolist(); props=client.compound_properties(cids,refresh) if cids else pd.DataFrame(); mapping=props.set_index("CID").to_dict("index") if not props.empty else {}; rows=[]
    for _,row in frame.iterrows():
        cid=int(row[cid_col]) if pd.notna(row.get(cid_col)) else None; prop=mapping.get(cid,{}) if cid else {}; smiles=row.get(smiles_col) or prop.get("SMILES") or prop.get("ConnectivitySMILES")
        try: structure=standardize_structure(smiles)
        except (StructureStandardizationError,TypeError): structure=None
        rows.append({"pubchem_cid":cid,"resolved_smiles":smiles,"pubchem_canonical_smiles":prop.get("ConnectivitySMILES") or prop.get("CanonicalSMILES"),"pubchem_isomeric_smiles":prop.get("SMILES") or prop.get("IsomericSMILES"),"pubchem_inchi":prop.get("InChI"),"pubchem_inchikey":prop.get("InChIKey"),"molecular_formula":prop.get("MolecularFormula"),"structure_valid":structure is not None,**(structure or {})})
    return pd.concat([frame.reset_index(drop=True),pd.DataFrame(rows)],axis=1)
