PROHIBITED=["successful drug discovered","new inhibitors were discovered","therapeutic candidate","clinical candidate"]
def claims_for_phases(completed):
 allowed=[]
 if 1 in completed:allowed.append("An explainable Ki-regression framework was developed and internally validated.")
 if 2 in completed:allowed.append("Prospective candidates were computationally prioritized.")
 if 3 in completed:allowed.append("Reproducible biochemical hits were identified from QC-passing laboratory evidence.")
 if 4 in completed:allowed.append("Orthogonally supported binders were identified from independent QC-passing evidence.")
 if 5 in completed:allowed.append("Cellular or early-lead evidence was classified according to imported QC-passing results.")
 return {"allowed":allowed,"prohibited":PROHIBITED}
