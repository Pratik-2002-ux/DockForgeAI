import pandas as pd,pytest
from src.features.fingerprint_features import fingerprint_matrix
def test_morgan_dimension(): assert fingerprint_matrix(pd.DataFrame({"canonical_isomeric_smiles":["CCO"]})).shape==(1,2048)
def test_invalid_fingerprint_molecule():
 with pytest.raises(ValueError): fingerprint_matrix(pd.DataFrame({"canonical_isomeric_smiles":["bad["]}))
