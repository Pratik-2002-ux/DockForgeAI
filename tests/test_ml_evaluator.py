import numpy as np
from src.ml.metrics import classification_metrics,expected_calibration_error
def test_metrics_include_primary_values(): assert {"roc_auc","pr_auc","mcc"}<classification_metrics(np.array([0,1]),np.array([.1,.9])).keys()
def test_ece_is_zero_for_perfect_confidence(): assert expected_calibration_error(np.array([0,1]),np.array([0.,1.]))==0
