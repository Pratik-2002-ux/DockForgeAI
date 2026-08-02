"""Scientific acceptance rules for PubChem HIV-1 protease assays."""
from __future__ import annotations

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def is_chembl_depositor(value): return "chembl" in str(value).lower()
def direct_protease_measurement(name,description,target=""):
    text=" ".join(map(str,(name,description,target))).lower(); has_target=("hiv-1 protease" in text or "hiv 1 protease" in text or ("gag-pol" in text and "protease" in text)); direct=any(x in text for x in ("enzyme inhibition","protease inhibition","inhibition of hiv","purified hiv","recombinant hiv","protease activity","inhibition constants")); disallowed=any(x in text for x in ("reverse transcriptase","integrase","cytotoxicity","pharmacokinetic")); return bool(has_target and direct and not disallowed)
def audit_assay(record,minimum_usable=30):
    aid=int(record.get("aid",record.get("uid",0))); depositor=record.get("currentsourcename") or ",".join(record.get("sourcenamelist",[])); name=record.get("assayname",""); description=record.get("assaydescription",""); targets="; ".join(x.get("name","") for x in record.get("proteintargetlist",[])); total=int(record.get("totalsidcount") or 0); active=int(record.get("activesidcount") or 0); inactive=max(total-active,0); reasons=[]
    independent=not is_chembl_depositor(depositor)
    direct=direct_protease_measurement(name,description,targets)
    if aid==162060 or not independent: reasons.append("depositor_is_chembl_or_chembl_derived")
    if not direct: reasons.append("not_clear_direct_hiv1_protease_measurement")
    if not (active>0 and inactive>0): reasons.append("active_and_inactive_outcomes_not_both_present")
    if total<minimum_usable: reasons.append(f"fewer_than_{minimum_usable}_compounds_before_overlap_removal")
    status="acceptable" if not reasons else "unsuitable"
    return {"pubchem_aid":aid,"assay_name":name,"depositor":depositor,"description":description,"target_name":targets,"assay_type":record.get("activityoutcomemethod",None),"total_tested_compounds":total,"active_count":active,"inactive_count":inactive,"inconclusive_count":None,"probe_count":record.get("probesidcount") or 0,"update_date":record.get("modifydate"),"depositor_independent_from_chembl":independent,"direct_target_measurement":direct,"experimental_outcomes_present":total>0,"active_and_inactive_present":active>0 and inactive>0,"result_definition_clear":direct,"structure_mapping_available":None,"training_exact_overlap_count":None,"training_scaffold_overlap_count":None,"usable_compound_count":total,"assay_quality_status":status,"rejection_reasons":"; ".join(reasons),"source_endpoint":f"{EUTILS}/esummary.fcgi?db=pcassay&id={aid}&retmode=json"}
