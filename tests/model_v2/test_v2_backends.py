import pytest
from src.models_v2.backends import availability,create_model,BackendUnavailable,BACKENDS
def test_all_requested_backends_registered(): assert len(BACKENDS)==10
def test_random_forest_available(): assert availability()["random_forest"]
@pytest.mark.parametrize("name",["chemprop_mpnn","gin","gine","attentivefp","graphgps","gin_transformer"])
def test_missing_graph_backend_is_explicit(name):
 if not availability()[name]:
  with pytest.raises(BackendUnavailable): create_model(name)
