"""Orchestrate all offline validation analyses."""
from __future__ import annotations
from .core import ensure_dirs
from .source_bias import run_source_bias
from .adversarial_validation import run_adversarial_validation
from .permutation_tests import run_permutations
from .dummy_baselines import run_dummy_baselines
from .nearest_neighbor_baseline import run_nearest_neighbor
from .similarity_stratification import run_similarity_stratification
from .scaffold_stratification import run_scaffold_stratification
from .descriptor_matching import run_descriptor_matching
from .reweighting import run_reweighting
from .threshold_sensitivity import run_threshold_sensitivity
from .decoy_sensitivity import run_decoy_sensitivity
from .calibration_robustness import run_calibration_robustness
from .stability_analysis import run_stability
from .explanation_sanity import run_explanation_sanity
from .probability_shift import run_probability_shift
from .validation_report import build_scorecard, generate_figures, generate_reports
from .validation_manifest import build_manifest

def run_all(permutations=100,bootstrap_iterations=500,decoy_seeds=range(2026,2036),report_only=False):
    import warnings
    from sklearn.exceptions import ConvergenceWarning
    warnings.filterwarnings("ignore", category=ConvergenceWarning)
    ensure_dirs()
    if not report_only:
        run_source_bias(); run_adversarial_validation(); run_permutations(permutations); run_dummy_baselines(); run_nearest_neighbor(); run_similarity_stratification(); run_scaffold_stratification(); run_descriptor_matching(); run_reweighting(); run_threshold_sensitivity(); run_decoy_sensitivity(decoy_seeds); run_calibration_robustness(bootstrap_iterations); run_stability(); run_explanation_sanity(); run_probability_shift()
    build_scorecard(); generate_figures(); report=generate_reports(); manifest=build_manifest(); return {"report":report,"artifact_count":len(manifest["artifacts"])}
