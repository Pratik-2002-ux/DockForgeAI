from __future__ import annotations
from scipy.stats import spearmanr
def importance_stability(a,b,top=20):
 da={r["feature"]:r["importance"] for r in a};db={r["feature"]:r["importance"] for r in b};common=sorted(set(da)&set(db)); correlation=float(spearmanr([da[k] for k in common],[db[k] for k in common]).statistic) if len(common)>1 else 0.; overlap=len({r["feature"] for r in a[:top]}&{r["feature"] for r in b[:top]})/top
 if abs(correlation-1.0)<1e-12: correlation=1.0
 return {"spearman_correlation":correlation,"top_feature_overlap":overlap}
