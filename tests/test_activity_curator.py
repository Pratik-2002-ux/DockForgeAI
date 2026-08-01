from src.datasets.activity_curator import curate_activities

CFG={"target_chembl_id":"CHEMBL243","standard_activity_types":["IC50"],"active_pchembl_threshold":6,"ambiguous_lower_threshold":5,"minimum_target_confidence_score":8}
def row(p, ident="1"): return {"activity_id":ident,"target_chembl_id":"CHEMBL243","standard_type":"IC50","standard_relation":"=","standard_value":str(10**(9-p)),"standard_units":"nM","pchembl_value":str(p),"canonical_smiles":"CCO","molecule_chembl_id":"CHEMBL1"}
def test_threshold_and_ambiguous_exclusion() -> None:
    result=curate_activities([row(6.0),row(5.5,"2")],CFG)
    assert len(result["active"])==1 and len(result["ambiguous"])==1
def test_conflicting_measurements_excluded() -> None:
    result=curate_activities([row(7),row(4,"2")],CFG)
    assert len(result["conflicting"])==1

def test_chembl_validity_flag_excluded() -> None:
    flagged=row(8); flagged["data_validity_comment"]="Outside typical range"
    assert curate_activities([flagged],CFG)["invalid"][0]["exclusion_reason"]=="outside_typical_range"
