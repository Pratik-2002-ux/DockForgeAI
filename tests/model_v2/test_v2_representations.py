import pytest
from src.models_v2.representations import fingerprints,graph_from_smiles
@pytest.mark.parametrize("kind,width",[("morgan",128),("ecfp6",128),("maccs",167)])
def test_fingerprint_dimensions(kind,width):
 x=fingerprints(["CCO"],kind,n_bits=128); assert x.shape==(1,width)
@pytest.mark.parametrize("smiles",["CCO","c1ccccc1","CC(=O)N","C[C@H](O)Cl"])
def test_graph_shapes(smiles):
 g=graph_from_smiles(smiles); assert g["atom_features"].shape[1]==5 and g["edge_index"].shape[0]==2 and g["bond_features"].shape[1]==4
