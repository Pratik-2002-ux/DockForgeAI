import numpy as np,pytest
from src.models_v2.robustness import perturb_descriptors,randomized_smiles,canonical_tautomer,remove_stereochemistry
def test_perturbation_deterministic(): assert np.array_equal(perturb_descriptors([[1]],seed=2),perturb_descriptors([[1]],seed=2))
@pytest.mark.parametrize("smiles",["CCO","CCN","c1ccccc1"])
def test_random_smiles_parse(smiles):
 from rdkit import Chem
 assert all(Chem.MolFromSmiles(x) for x in randomized_smiles(smiles,3))
def test_tautomer_returns_smiles(): assert canonical_tautomer("CC(=O)C").strip()
def test_stereo_removed(): assert "@" not in remove_stereochemistry("C[C@H](O)Cl")
