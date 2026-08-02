import pytest
from src.models_v2.curation import to_nm,p_activity,assay_quality_score
@pytest.mark.parametrize("value,unit,expected",[(1,"nM",1),(1,"uM",1000),(1,"mM",1e6),(1,"pM",.001)])
def test_unit_conversion(value,unit,expected): assert to_nm(value,unit)==expected
@pytest.mark.parametrize("nm,p",[(1,9),(10,8),(100,7),(1000,6)])
def test_p_activity(nm,p): assert p_activity(nm)==p
@pytest.mark.parametrize("confidence",[0,2,5,8,10])
def test_quality_is_bounded(confidence): assert 0<=assay_quality_score(confidence)<=1
def test_unsupported_unit_rejected():
 with pytest.raises(ValueError): to_nm(1,"mg")
