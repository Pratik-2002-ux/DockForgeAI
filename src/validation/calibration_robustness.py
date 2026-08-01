"""Calibration robustness across domain and chemical strata."""
from __future__ import annotations
import numpy as np, pandas as pd
from scipy.special import logit
from sklearn.linear_model import LogisticRegression
from .core import OUTPUT, best_predictions, safe_metrics, write_json

def calibration_slope_intercept(y,p):
    z=logit(np.clip(np.asarray(p),1e-6,1-1e-6)).reshape(-1,1); model=LogisticRegression(C=1e6,max_iter=1000).fit(z,y); return float(model.coef_[0,0]),float(model.intercept_[0])
def maximum_calibration_error(y,p,bins=10):
    frame=pd.DataFrame({"y":y,"p":p,"bin":pd.cut(p,np.linspace(0,1,bins+1),include_lowest=True)}); gaps=[abs(g.y.mean()-g.p.mean()) for _,g in frame.groupby("bin",observed=True)]; return float(max(gaps,default=0))
def run_calibration_robustness(bootstrap_iterations=500, output=OUTPUT/"calibration"):
    pred=best_predictions(); validation=best_predictions("validation"); y=pred.true_label.eq("active").astype(int).to_numpy(); p=pred.predicted_probability.to_numpy(); vy=validation.true_label.eq("active").astype(int); calibrator=LogisticRegression(C=1e6,max_iter=1000).fit(logit(np.clip(validation.raw_probability,1e-6,1-1e-6)).to_numpy().reshape(-1,1),vy); platt=calibrator.predict_proba(logit(np.clip(pred.raw_probability,1e-6,1-1e-6)).to_numpy().reshape(-1,1))[:,1]; ad=pd.read_csv(OUTPUT.parent/"xai_predictions/applicability_domain.csv"); pred=pred.merge(ad[["compound_id","in_domain"]],on="compound_id"); pred["raw"]=pred.raw_probability; pred["platt"]=platt; pred["isotonic"]=p; rows=[]; slope=[]
    for method in ("raw","platt","isotonic"):
      for group,frame in [("all",pred),("in_domain",pred[pred.in_domain]),("out_of_domain",pred[~pred.in_domain])]:
        yy=frame.true_label.eq("active").astype(int); pp=frame[method]; m=safe_metrics(yy,pp); m.update(group=group,calibration=method,maximum_calibration_error=maximum_calibration_error(yy,pp)); rows.append(m)
        if yy.nunique()==2: s,i=calibration_slope_intercept(yy,pp); slope.append({"group":group,"calibration":method,"slope":s,"intercept":i})
    rng=np.random.default_rng(2026); briers=[]
    for _ in range(bootstrap_iterations):
        ids=rng.integers(0,len(y),len(y)); briers.append(np.mean((y[ids]-p[ids])**2))
    pd.DataFrame(rows).to_csv(output/"calibration_robustness.csv",index=False); pd.DataFrame(slope).to_csv(output/"calibration_slope_intercept.csv",index=False)
    bin_rows=[]
    for method in ("raw","platt","isotonic"):
        b=pd.DataFrame({"true_label":y,"probability":pred[method],"bin":pd.cut(pred[method],np.linspace(0,1,11),include_lowest=True).astype(str)}).groupby("bin",observed=True).agg(count=("true_label","size"),observed_rate=("true_label","mean"),mean_probability=("probability","mean")).reset_index(); b["calibration"]=method; bin_rows.append(b)
    pd.concat(bin_rows).to_csv(output/"calibration_bin_data.csv",index=False)
    summary={"bootstrap_iterations":bootstrap_iterations,"brier_ci_95":[float(np.quantile(briers,.025)),float(np.quantile(briers,.975))],"calibration_selected_on_validation_only":True,"groups":rows}; write_json(output/"calibration_summary.json",summary); return pd.DataFrame(rows)
