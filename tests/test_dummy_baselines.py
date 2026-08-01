import numpy as np
from src.validation.core import safe_metrics
def test_constant_baseline_is_random_auc(): assert safe_metrics([0,1],[.5,.5])["roc_auc"]==.5
def test_metric_confusion_counts_complete():
    m=safe_metrics(np.array([0,1]),np.array([.2,.8])); assert sum(m[k] for k in ("tn","fp","fn","tp"))==2
