#!/usr/bin/env python
"""Run the offline XAI validation suite."""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))
from src.validation.runner import run_all

def parser():
    p=argparse.ArgumentParser(); modes=p.add_mutually_exclusive_group(); modes.add_argument("--all",action="store_true"); modes.add_argument("--source-bias-only",action="store_true"); modes.add_argument("--permutation-only",action="store_true"); modes.add_argument("--baseline-only",action="store_true"); modes.add_argument("--matched-subset-only",action="store_true"); modes.add_argument("--decoy-sensitivity-only",action="store_true"); modes.add_argument("--calibration-only",action="store_true"); modes.add_argument("--explanation-sanity-only",action="store_true"); modes.add_argument("--report-only",action="store_true"); p.add_argument("--overwrite",action="store_true"); p.add_argument("--strict",action="store_true"); p.add_argument("--offline",action="store_true"); p.add_argument("--permutations",type=int,default=100); p.add_argument("--bootstrap-iterations",type=int,default=500); p.add_argument("--decoy-seeds",type=int,nargs="+",default=list(range(2026,2036))); return p
def main():
    a=parser().parse_args()
    if not a.offline: raise SystemExit("Validation must be run with --offline")
    if a.all or a.report_only or not any((a.source_bias_only,a.permutation_only,a.baseline_only,a.matched_subset_only,a.decoy_sensitivity_only,a.calibration_only,a.explanation_sanity_only)):
        result=run_all(a.permutations,a.bootstrap_iterations,a.decoy_seeds,a.report_only)
    else:
        from src.validation.core import ensure_dirs; ensure_dirs(); result={"mode":"selected","offline":True}
        if a.source_bias_only:
            from src.validation.source_bias import run_source_bias; from src.validation.adversarial_validation import run_adversarial_validation; run_source_bias(); run_adversarial_validation()
        elif a.permutation_only:
            from src.validation.permutation_tests import run_permutations; run_permutations(a.permutations)
        elif a.baseline_only:
            from src.validation.dummy_baselines import run_dummy_baselines; from src.validation.nearest_neighbor_baseline import run_nearest_neighbor; run_dummy_baselines(); run_nearest_neighbor()
        elif a.matched_subset_only:
            from src.validation.descriptor_matching import run_descriptor_matching; run_descriptor_matching()
        elif a.decoy_sensitivity_only:
            from src.validation.decoy_sensitivity import run_decoy_sensitivity; run_decoy_sensitivity(a.decoy_seeds)
        elif a.calibration_only:
            from src.validation.calibration_robustness import run_calibration_robustness; run_calibration_robustness(a.bootstrap_iterations)
        elif a.explanation_sanity_only:
            from src.validation.explanation_sanity import run_explanation_sanity; run_explanation_sanity()
    print(json.dumps(result,indent=2))
if __name__=="__main__": main()
