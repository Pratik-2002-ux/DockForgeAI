from src.validation.descriptor_matching import standardized_mean_difference, MATCH_COLUMNS
def test_equal_groups_have_zero_smd(): assert standardized_mean_difference([1,2],[1,2])==0
def test_matching_uses_required_descriptors(): assert set(["molecular_weight","logp","tpsa","hbd","hba","rotatable_bonds","formal_charge","heavy_atom_count"])==set(MATCH_COLUMNS)
