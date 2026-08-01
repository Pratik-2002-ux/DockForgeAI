import pytest
from src.features.feature_schema import validate_feature_names
def test_safe_feature_names(): validate_feature_names(["logp","Morgan_1"])
def test_leakage_feature_rejected():
 with pytest.raises(ValueError): validate_feature_names(["source_database"])
