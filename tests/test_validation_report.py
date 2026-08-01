from src.validation.validation_report import LIMITATIONS
def test_report_states_tests_are_not_biological_truth(): assert any("not biological truth" in x for x in LIMITATIONS)
def test_report_states_decoy_limitation(): assert any("not experimentally verified inactive" in x for x in LIMITATIONS)
