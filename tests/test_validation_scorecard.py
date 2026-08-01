from src.validation.validation_report import LIMITATIONS
def test_external_validation_not_implied(): assert any("has not yet been completed" in x for x in LIMITATIONS)
def test_high_performance_not_therapeutic_claim(): assert any("does not prove therapeutic" in x for x in LIMITATIONS)
