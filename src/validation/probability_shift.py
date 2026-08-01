"""Probability distribution shift diagnostics."""
from __future__ import annotations
import numpy as np, pandas as pd
from scipy.spatial.distance import jensenshannon
from scipy.stats import ks_2samp, wasserstein_distance
from .core import OUTPUT, best_predictions, load_data

def run_probability_shift(output=OUTPUT/"probability_shift"):
    d=load_data(); p=best_predictions(); frame=d.merge(p[["compound_id","predicted_probability","predicted_class"]],on="compound_id"); summaries=[]; distances=[]
    for field in ("source_database","label_source","dataset_split","binary_label"):
        groups=[]
        for name,g in frame.groupby(field): summaries.append({"field":field,"group":str(name),"count":len(g),"mean_probability":g.predicted_probability.mean(),"std_probability":g.predicted_probability.std()}); groups.append((str(name),g.predicted_probability.to_numpy()))
        for i,(a,x) in enumerate(groups):
            for b,y in groups[i+1:]:
                hx,_=np.histogram(x,bins=20,range=(0,1),density=True); hy,_=np.histogram(y,bins=20,range=(0,1),density=True); hx=(hx+1e-12)/sum(hx+1e-12); hy=(hy+1e-12)/sum(hy+1e-12)
                distances.append({"field":field,"group_a":a,"group_b":b,"ks_statistic":ks_2samp(x,y).statistic,"wasserstein_distance":wasserstein_distance(x,y),"jensen_shannon_divergence":jensenshannon(hx,hy)**2})
    pd.DataFrame(summaries).to_csv(output/"probability_shift_summary.csv",index=False); pd.DataFrame(distances).to_csv(output/"distribution_distances.csv",index=False); return pd.DataFrame(distances)
