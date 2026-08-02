from src.external_validation.pubchem_external_curator import OUTCOME_MAP
def test_only_supported_outcomes_are_binary(): assert set(OUTCOME_MAP)=={"active","inactive"}
def test_no_threshold_is_defined_in_outcome_map(): assert all(isinstance(x,str) for x in OUTCOME_MAP.values())
