ALLOWED=["computational candidate","biochemical hit","orthogonally supported binder","cellularly active hit","selective cellular hit","early lead","lead with unresolved liabilities","not advanced","inconclusive"]
def maximum_status(evidence):
 if evidence.get("early_safety_qc") and evidence.get("adme_qc") and evidence.get("selective_cellular"):return "early lead"
 if evidence.get("selective_cellular"):return "selective cellular hit"
 if evidence.get("cellular_qc"):return "cellularly active hit"
 if evidence.get("orthogonal_qc"):return "orthogonally supported binder"
 if evidence.get("biochemical_qc"):return "biochemical hit"
 return "computational candidate"
