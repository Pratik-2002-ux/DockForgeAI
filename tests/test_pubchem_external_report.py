from src.external_validation.pubchem_external_report import STATEMENTS,SECTIONS
def test_report_states_locked_model(): assert any("remained locked" in x for x in STATEMENTS)
def test_report_rejects_aid_162060(): assert any("162060" in x and "ChEMBL" in x for x in STATEMENTS)
