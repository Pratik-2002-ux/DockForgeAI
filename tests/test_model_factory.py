import pytest
from src.ml.model_factory import create_model
def test_all_required_models_construct(): assert all(create_model(n,2026) is not None for n in ("logistic_regression","random_forest","extra_trees","hist_gradient_boosting"))
def test_unknown_model_rejected():
 with pytest.raises(ValueError): create_model("unknown",1)
