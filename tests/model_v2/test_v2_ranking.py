import pandas as pd
from src.models_v2.ranking import rank_models,DEFAULT_WEIGHTS
def data(): return pd.DataFrame({"model":["a","b"],"roc_auc":[.9,.8],"pr_auc":[.9,.8],"mcc":[.8,.7],"balanced_accuracy":[.8,.7],"calibration":[.1,.2],"external_validation":[.8,.7],"novel_scaffold":[.8,.7],"inference_speed":[1,2],"memory_usage":[1,2]})
def test_weights_sum_one(): assert abs(sum(DEFAULT_WEIGHTS.values())-1)<1e-12
def test_better_model_ranks_first(): assert rank_models(data()).iloc[0].model=="a"
def test_input_not_mutated():
 d=data(); rank_models(d); assert "selection_score" not in d
