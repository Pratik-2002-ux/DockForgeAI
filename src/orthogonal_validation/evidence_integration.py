def integrate(biochemical_status,orthogonal_qc,orthogonal_positive):
 if biochemical_status not in {"reproducible biochemical hit","quantitative biochemical hit"}:return "no orthogonal evidence"
 if orthogonal_qc!="pass":return "orthogonal test failed"
 return "orthogonally supported binder" if orthogonal_positive else "conflicting evidence"
