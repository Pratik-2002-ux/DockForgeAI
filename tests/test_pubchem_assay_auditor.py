from src.external_validation.pubchem_assay_auditor import audit_assay
def record(aid=1,source="BindingDB",total=40,active=20): return {"aid":aid,"currentsourcename":source,"totalsidcount":total,"activesidcount":active,"assayname":"HIV-1 protease inhibition","assaydescription":"Recombinant HIV-1 protease activity enzyme inhibition"}
def test_chembl_depositor_rejected(): assert audit_assay(record(source="ChEMBL"))["assay_quality_status"]=="unsuitable"
def test_aid_162060_always_rejected(): assert "chembl" in audit_assay(record(aid=162060))["rejection_reasons"]
