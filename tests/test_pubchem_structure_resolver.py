import pandas as pd
from src.external_validation.pubchem_structure_resolver import resolve_result_structures
class Client:
 def compound_properties(self,*a,**k): return pd.DataFrame([{"CID":1,"SMILES":"CCO","InChI":"x","InChIKey":"y","MolecularFormula":"C2H6O"}])
def test_sid_cid_structure_resolution(): assert resolve_result_structures(Client(),pd.DataFrame([{"PUBCHEM_CID":1,"PUBCHEM_EXT_DATASOURCE_SMILES":"CCO"}])).structure_valid.item()
def test_invalid_structure_flagged(): assert not resolve_result_structures(Client(),pd.DataFrame([{"PUBCHEM_CID":2,"PUBCHEM_EXT_DATASOURCE_SMILES":"invalid"}])).structure_valid.item()
