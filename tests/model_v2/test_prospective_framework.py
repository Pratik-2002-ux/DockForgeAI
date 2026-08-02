import math
import pytest
from src.models_v2.regression.ki_parser import parse_ki,pki_from_nm
from src.models_v2.regression.censored_measurements import is_censored
from src.models_v2.regression.regression_models import models
from src.experimental_validation.hit_classification import classify
from src.orthogonal_validation.evidence_integration import integrate
from src.orthogonal_validation.study_schema import validate_method
from src.lead_validation.evidence_grading import maximum_status
from src.lead_validation.result_ingestion import selectivity_index
from src.prospective_screening.provenance import immutable_candidate_id
from src.publication_claims import claims_for_phases

@pytest.mark.parametrize("value,relation,nm",[("1","=",1),("<2","<",2),("<=3","<=",3),(">4",">",4),(">=5",">=",5),("1e2","=",100),("0",None,None),("bad",None,None)])
def test_ki_parsing(value,relation,nm):assert parse_ki(value)==(relation,nm)
@pytest.mark.parametrize("nm,expected",[(1,9),(10,8),(100,7),(1000,6),(10000,5),(1e6,3)])
def test_pki_conversion(nm,expected):assert pki_from_nm(nm)==pytest.approx(expected)
@pytest.mark.parametrize("relation,expected",[("=",False),("<",True),("<=",True),(">",True),(">=",True)])
def test_censored_detection(relation,expected):assert is_censored(relation)==expected
@pytest.mark.parametrize("qc,rep,quant,signal,expected",[(False,False,False,False,"test failed"),(True,True,False,False,"no confirmed inhibition"),(True,False,False,True,"inconclusive"),(True,True,False,True,"reproducible biochemical hit"),(True,True,True,True,"quantitative biochemical hit"),(False,True,True,True,"test failed"),(True,False,True,True,"inconclusive"),(True,True,True,False,"no confirmed inhibition")])
def test_biochemical_status(qc,rep,quant,signal,expected):assert classify(qc,rep,quant,signal)==expected
@pytest.mark.parametrize("bio,qc,pos,expected",[("no confirmed inhibition","pass",True,"no orthogonal evidence"),("reproducible biochemical hit","fail",True,"orthogonal test failed"),("reproducible biochemical hit","pass",True,"orthogonally supported binder"),("reproducible biochemical hit","pass",False,"conflicting evidence"),("quantitative biochemical hit","pass",True,"orthogonally supported binder"),("inconclusive","pass",True,"no orthogonal evidence")])
def test_orthogonal_integration(bio,qc,pos,expected):assert integrate(bio,qc,pos)==expected
@pytest.mark.parametrize("evidence,expected",[({},"computational candidate"),({"biochemical_qc":True},"biochemical hit"),({"orthogonal_qc":True},"orthogonally supported binder"),({"cellular_qc":True},"cellularly active hit"),({"selective_cellular":True},"selective cellular hit"),({"selective_cellular":True,"adme_qc":True,"early_safety_qc":True},"early lead")])
def test_evidence_status(evidence,expected):assert maximum_status(evidence)==expected
@pytest.mark.parametrize("eff,tox,expected",[(1,10,10),(2,10,5),(5,5,1),(.5,5,10),(10,100,10),(4,20,5),(3,9,3),(8,16,2),(.25,1,4),(7,21,3)])
def test_selectivity_index(eff,tox,expected):assert selectivity_index(eff,tox,"nM","nM")==pytest.approx(expected)
@pytest.mark.parametrize("source",[f"C{i}" for i in range(10)])
def test_candidate_id_stable(source):assert immutable_candidate_id("STUDY",source)==f"STUDY:{source}"
@pytest.mark.parametrize("phases",[{1},{1,2},{1,2,3},{1,2,3,4},{1,2,3,4,5}])
def test_publication_claims_never_claim_drug(phases):assert "successful drug discovered" in claims_for_phases(phases)["prohibited"]
@pytest.mark.parametrize("method",["surface_binding","calorimetric_binding","thermophoretic_binding","thermal_shift","other_justified"])
def test_orthogonal_methods(method):assert validate_method(method)==method
@pytest.mark.parametrize("name",["dummy","linear","ridge","elastic_net","random_forest","extra_trees","hist_gradient_boosting"])
def test_mandatory_regressors(name):assert name in models(2026)
