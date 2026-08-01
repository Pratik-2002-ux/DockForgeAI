from src.validation.core import safe_metrics
def test_single_class_stratum_auc_is_undefined(): assert safe_metrics([1,1],[.6,.8])["roc_auc"] is None
def test_similarity_metric_counts_rows(): assert safe_metrics([0,1],[.1,.9])["sample_count"]==2
