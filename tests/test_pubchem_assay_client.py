import pytest
from src.external_validation.pubchem_assay_client import PubChemAssayClient,PubChemRequestError
def test_offline_cache_mode_rejects_miss(tmp_path):
    with pytest.raises(PubChemRequestError): PubChemAssayClient(tmp_path,offline=True).get("https://example.invalid")
def test_cached_response_is_reused(tmp_path):
    c=PubChemAssayClient(tmp_path,offline=True); c._path("u").write_bytes(b"x"); assert c.get("u")==b"x"
