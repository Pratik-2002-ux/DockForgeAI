from src.external_validation.pubchem_external_evaluator import group_metrics
import pandas as pd
def test_sparse_group_marked_unstable(): assert group_metrics(pd.DataFrame({"true_label":["active","inactive_or_decoy"],"predicted_probability":[.8,.2]}),10)["unstable"]
def test_group_metrics_preserve_counts(): assert group_metrics(pd.DataFrame({"true_label":["active","inactive_or_decoy"],"predicted_probability":[.8,.2]}),10)["sample_count"]==2
