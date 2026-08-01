from src.validation.permutation_tests import empirical_p_value
def test_empirical_p_value_has_plus_one_correction(): assert empirical_p_value([.1,.2,.3],.4)==.25
def test_empirical_p_value_counts_ties(): assert empirical_p_value([.5,.6],.5)==1.0
