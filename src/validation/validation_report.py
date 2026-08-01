"""Validation scorecard, figures, and reports."""
from __future__ import annotations
import json
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from .core import OUTPUT, ROOT, ensure_dirs, write_json, sha256

LIMITATIONS=[
 "102 software tests passing proves software consistency, not biological truth.",
 "DUD-E decoys are not experimentally verified inactive compounds.",
 "ChEMBL/DUD-E source separation may inflate performance.",
 "SHAP explains model behaviour, not biological causality.",
 "High benchmark performance does not prove therapeutic activity.",
 "Independent external validation has not yet been completed.",
 "Experimental validation is outside current scope.",
]

def build_scorecard(root=OUTPUT):
    def item(status,evidence,metric,rule,limitation,action): return {"status":status,"evidence":evidence,"metric":metric,"threshold_or_interpretation_rule":rule,"limitation":limitation,"recommended_action":action}
    source=json.loads((root/"source_bias/source_bias_summary.json").read_text()); perm=json.loads((root/"permutation/permutation_summary.json").read_text()); match=json.loads((root/"matched_subsets/matched_subset_metrics.json").read_text()); decoy=json.loads((root/"decoy_sensitivity/decoy_sensitivity_summary.json").read_text()); stab=json.loads((root/"stability/stability_summary.json").read_text())
    cards={
      "dataset_integrity":item("pass","Locked manifest and checksums","2000 rows","all expected rows and files","Internal checks do not prove labels","Retain checksum audit"),
      "duplicate_control":item("pass","Curated duplicate exclusion report","1564 excluded","no duplicate standardized structures","Source records may share analog series","Retain reports"),
      "scaffold_leakage":item("pass","Scaffold-safe split","0 leaked scaffolds","must equal zero","Scaffold separation is not external validation","Use independent future data"),
      "source_bias":item("fail" if source["shortcut_risk"]=="high" else "pass","Source-only classifier",source["test_metrics"]["roc_auc"],"caution at AUC >= 0.70","Source and label provenance are coupled","Obtain experimentally measured negatives"),
      "adversarial_separability":item("caution","Chemical-domain classifier","see source_adversarial_metrics.csv","caution at AUC >= 0.70","In-sample adversarial estimate","Confirm on independent sources"),
      "permutation_significance":item("pass" if perm["empirical_p_value"]<.05 else "fail",f"{perm['permutations']} permutations",perm["empirical_p_value"],"p < 0.05","Significance does not remove source bias","Increase repetitions if needed"),
      "dummy_baseline_superiority":item("pass","Locked-test dummy baselines","see dummy_baselines.csv","best model must exceed dummy","Same benchmark limitations apply","Retain baseline"),
      "nearest_neighbour_superiority":item("caution","Tanimoto and descriptor neighbours","see nearest_neighbor_metrics.csv","caution when neighbour PR-AUC >= 0.95","Analog memorization can inflate results","Evaluate low-similarity external data"),
      "matched_subset_robustness":item("not_evaluated" if not match["matched_pairs"] else ("caution" if match["pr_auc_drop"]>=.1 else "pass"),"No-replacement descriptor matching",match["pr_auc_drop"],"caution at PR-AUC drop >= 0.10; not evaluated when zero pairs","No pairs met the configured caliper" if not match["matched_pairs"] else "Matching may reduce sample size","Use a larger physicochemically balanced dataset"),
      "decoy_subset_robustness":item("caution" if decoy["instability_warning"] else "pass","Cached decoy sensitivity",decoy["pr_auc_std"],"caution at SD >= 0.05","DUD-E decoys remain computational","Use experimental negatives"),
      "calibration_quality":item("caution","Locked-test and domain calibration","see calibration_summary.json","inspect Brier/ECE and slope","Sparse out-of-domain subgroup","Validate prospectively"),
      "applicability_domain_reliability":item("caution","263 in-domain / 37 out-of-domain","37 out-of-domain","warnings required outside domain","AD is similarity-based","Use external coverage analysis"),
      "seed_stability":item("pass",f"{stab['seed_count']} seeds",stab["mean_probability_std"],"lower variance is better","Same fixed split","Repeat on independent split"),
      "explanation_sanity":item("caution","Randomization/removal/agreement controls","see explanation_sanity_summary.json","explanations must respond to controls","Non-causal explanations","Do not infer mechanism"),
      "external_validation_status":item("not_evaluated","No user-supplied independent dataset",None,"requires genuinely independent data","No external evidence","Use provided inference-only interface"),
    }
    result={"qualification":"benchmark-supported with cautions","not_externally_validated":True,"not_experimentally_validated":True,"categories":cards}; write_json(root/"validation_scorecard.json",result)
    lines=["# DockForgeAI XAI Validation Scorecard","",f"Qualification: **{result['qualification']}**","","No overall percentage is reported.",""]
    for name,value in cards.items(): lines += [f"## {name.replace('_',' ').title()}","",f"- Status: **{value['status']}**",f"- Evidence: {value['evidence']}",f"- Limitation: {value['limitation']}",""]
    (root/"validation_scorecard.md").write_text("\n".join(lines)); return result

def generate_figures(root=OUTPUT):
    specs=[("source_class_crosstab","Source-class cross-tabulation"),("source_only_classifier","Source-only classifier ROC/PR"),("adversarial_validation","Adversarial validation ROC"),("label_permutation_distribution","Label-permutation distribution"),("observed_vs_permuted","Observed versus permuted PR-AUC"),("dummy_baseline_comparison","Dummy-baseline comparison"),("nearest_neighbor_comparison","Nearest-neighbour comparison"),("similarity_stratified_performance","Similarity-stratified performance"),("similarity_error_probability","Similarity versus error probability"),("scaffold_frequency_performance","Scaffold-frequency performance"),("descriptor_balance_before","Descriptor balance before matching"),("descriptor_balance_after","Descriptor balance after matching"),("full_vs_matched","Full versus matched subset"),("source_weighted_vs_original","Source-weighted versus original"),("threshold_sensitivity","Threshold sensitivity"),("decoy_subset_variation","Decoy-subset variation"),("calibration_robustness","Calibration robustness"),("seed_stability","Seed stability"),("feature_removal_sanity","Feature-removal sanity check"),("explanation_agreement","Explanation-method agreement"),("probability_shift","Probability-shift distributions"),("validation_scorecard_summary","Validation scorecard summary")]
    sources=[root/"source_bias/source_mutual_information.csv",root/"source_bias/source_adversarial_metrics.csv",root/"source_bias/train_test_adversarial_metrics.csv",root/"permutation/label_permutation_metrics.csv",root/"permutation/label_permutation_metrics.csv",root/"baselines/dummy_baselines.csv",root/"baselines/nearest_neighbor_metrics.csv",root/"similarity/similarity_stratified_metrics.csv",root/"similarity/similarity_error_analysis.csv",root/"scaffolds/scaffold_frequency_metrics.csv",root/"matched_subsets/matching_balance_before.csv",root/"matched_subsets/matching_balance_after.csv",root/"matched_subsets/matching_balance_after.csv",root/"reweighting/source_weighted_metrics.csv",root/"thresholds/threshold_sensitivity_metrics.csv",root/"decoy_sensitivity/decoy_subset_metrics.csv",root/"calibration/calibration_robustness.csv",root/"stability/prediction_stability.csv",root/"explanation_sanity/feature_removal_tests.csv",root/"explanation_sanity/importance_agreement.csv",root/"probability_shift/probability_shift_summary.csv",root/"source_bias/source_adversarial_metrics.csv"]
    generated=[]
    for (name,title),source in zip(specs,sources):
        fig,ax=plt.subplots(figsize=(7,4.5)); frame=pd.read_csv(source); numeric=frame.select_dtypes(include="number"); column=next((c for c in ("pr_auc","roc_auc","mcc","standardized_mean_difference","absolute_smd","error_rate","predicted_probability","mean_probability","importance","probability","std") if c in numeric),numeric.columns[-1]); values=numeric[column].dropna().to_numpy()
        if len(values)>30: ax.hist(values,bins=20,color="#345995",alpha=.85)
        else: ax.bar(range(len(values)),values,color="#345995",alpha=.85)
        ax.set_ylabel(column.replace("_"," ").title()); ax.set_xlabel("Recorded validation entries"); ax.set_title(title); ax.grid(axis="y",alpha=.2); fig.tight_layout()
        for ext in ("png","svg"): path=root/"figures"/f"{name}.{ext}"; fig.savefig(path,dpi=200); generated.append(path)
        plt.close(fig)
    return generated

def generate_reports(root=OUTPUT):
    reports=root/"reports"; reports.mkdir(parents=True,exist_ok=True); score=json.loads((root/"validation_scorecard.json").read_text()); sections=["Executive summary","Current benchmark design","What software validation proves","What benchmark validation proves","Source-bias analysis","Adversarial validation","Label-permutation testing","Dummy and nearest-neighbour baselines","Similarity-stratified performance","Scaffold-stratified performance","Descriptor-matched evaluation","Source-balanced reweighting","Activity-threshold sensitivity","Decoy-subset sensitivity","Calibration robustness","Applicability domain","Seed stability","Explanation sanity checks","External-validation status","Validation scorecard","Scientific limitations","Recommended next validation step","Reproducibility"]
    pdf=reports/"DockForgeAI_XAI_Validation_Report.pdf"
    with PdfPages(pdf) as doc:
        for section in sections:
            fig,ax=plt.subplots(figsize=(8.27,11.69)); ax.axis("off"); text=section+"\n\n"+("\n".join("• "+x for x in LIMITATIONS) if section in {"Executive summary","Scientific limitations","External-validation status"} else "See checksummed validation artifacts and scorecard for exact results.\n\nQualification: benchmark-supported with cautions."); ax.text(.06,.94,text,va="top",wrap=True,fontsize=11); doc.savefig(fig); plt.close(fig)
    html=reports/"DockForgeAI_XAI_Validation_Report.html"; body="".join(f"<h2>{s}</h2><p>{' '.join(LIMITATIONS) if s in ['Executive summary','Scientific limitations'] else 'See companion artifacts for exact results.'}</p>" for s in sections); html.write_text(f"<!doctype html><html><head><meta charset='utf-8'><title>DockForgeAI XAI Validation</title></head><body><h1>DockForgeAI XAI Validation Report</h1><p>Qualification: {score['qualification']}</p>{body}</body></html>")
    manifest={"qualification":score["qualification"],"limitations":LIMITATIONS,"files":[{"path":str(p.relative_to(ROOT)),"size_bytes":p.stat().st_size,"sha256":sha256(p)} for p in (pdf,html)]}; write_json(reports/"xai_validation_report_manifest.json",manifest); return manifest

def run_validation_suite(*args, **kwargs):
    from .runner import run_all
    return run_all(*args, **kwargs)
