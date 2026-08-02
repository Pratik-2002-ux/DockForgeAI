from src.external_validation import SEARCH_TERMS
def test_required_search_terms_present(): assert "HIV-1 protease" in SEARCH_TERMS and "retroviral protease inhibition" in SEARCH_TERMS
def test_search_terms_are_unique(): assert len(SEARCH_TERMS)==len(set(SEARCH_TERMS))
