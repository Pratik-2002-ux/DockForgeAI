import pandas as pd
import pytest

from src.models_v2.data.activity_harmonizer import classify_measurement, convert_to_nm, endpoint_fields
from src.models_v2.data.dataset_report import readiness
from src.models_v2.data.external_overlap import mark_external_overlap
from src.models_v2.data.measurement_parser import parse_measurement
from src.models_v2.data.pubchem_loader import load_pubchem
from src.models_v2.data.split_builder import add_splits, leakage_count


@pytest.mark.parametrize("text,relation,value",[("10","=",10),("<5","<",5),("<=7","<=",7),(">10000",">",10000),(">=2",">=",2),("1e3","=",1000)])
def test_relation_parsing(text,relation,value): assert parse_measurement(text)==(relation,value)

@pytest.mark.parametrize("value,units,expected",[(1,"nM",1),(1,"uM",1000),(1,"µM",1000),(1,"mM",1_000_000),(1,"pM",.001)])
def test_multisource_unit_conversion(value,units,expected): assert convert_to_nm(value,units)==expected

@pytest.mark.parametrize("value,relation,expected",[(1000,"=","active"),(1000,"<","active"),(1001,"=","ambiguous"),(9999,"=","ambiguous"),(10000,"=","inactive"),(10000,">=","inactive"),(None,"=","ambiguous")])
def test_binary_policy(value,relation,expected): assert classify_measurement(value,relation,endpoint="Ki")==expected

@pytest.mark.parametrize("outcome",["active","inactive"])
def test_explicit_outcomes(outcome): assert classify_measurement(None,None,outcome)==outcome

@pytest.mark.parametrize("endpoint,key",[("KI","ki_nm"),("KD","kd_nm"),("IC50","ic50_nm"),("EC50","ec50_nm")])
def test_endpoint_specific_fields(endpoint,key): assert endpoint_fields(endpoint,10)[key]==10

def test_unsupported_endpoint_is_ambiguous(): assert classify_measurement(1,"=",endpoint="kon")=="ambiguous"

def test_locked_pubchem_rejected(tmp_path):
    locked=tmp_path/"locked.csv"; locked.write_text("x\n")
    with pytest.raises(ValueError,match="locked"): load_pubchem(locked,locked)

def test_overlap_connectivity_flag():
    frame=pd.DataFrame({"canonical_isomeric_smiles":["CC"],"inchikey":["A"],"inchikey_connectivity":["BLOCK"],"murcko_scaffold_smiles":[""]})
    locked={"canonical_isomeric_smiles":set(),"inchikey":set(),"inchikey_connectivity":{"BLOCK"},"murcko_scaffold_smiles":set()}
    assert bool(mark_external_overlap(frame,locked).external_overlap_flag.iloc[0])

def test_all_group_splits_do_not_leak():
    frame=pd.DataFrame({"murcko_scaffold_smiles":[f"s{i//2}" for i in range(40)],"source_databases":[f"src{i%2}" for i in range(40)],"source_assay_ids":[f"a{i//3}" for i in range(40)],"publication_date":range(1980,2020)})
    out=add_splits(frame)
    assert leakage_count(out,"murcko_scaffold_smiles","dataset_split_scaffold")==0
    assert leakage_count(out,"source_databases","dataset_split_source")==0
    assert leakage_count(out,"source_assay_ids","dataset_split_assay")==0

def test_readiness_thresholds():
    base={"final_active_count":200,"final_inactive_count":200,"source_count":2,"assay_count":10,"unique_scaffold_count":100,"unique_compound_count":500,"locked_external_overlap_count":0}
    assert readiness(base)==(True,True)
    base["final_inactive_count"]=199
    assert readiness(base)==(False,False)
