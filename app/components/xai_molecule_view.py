from rdkit import Chem
from src.visualization.molecule_renderer import render_molecule_bytes
def molecule_png(smiles:str,compound_id:str)->bytes:
 mol=Chem.MolFromSmiles(smiles)
 if mol is None: raise ValueError("Invalid molecule")
 return render_molecule_bytes(mol,compound_id)
