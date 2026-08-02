from src.external_validation.pubchem_external_curator import map_outcome
def test_active_inactive_mapping(): assert map_outcome("Active")=="active" and map_outcome("Inactive")=="inactive_or_decoy"
def test_inconclusive_is_not_mapped(): assert map_outcome("Inconclusive") is None
