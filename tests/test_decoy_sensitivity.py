from src.validation.decoy_sensitivity import _properties
def test_cached_decoy_descriptor_shape(): assert len(_properties("CCO")[0])==8
def test_invalid_decoy_is_rejected(): assert _properties("not_smiles") is None
