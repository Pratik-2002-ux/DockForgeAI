from src.explainability.explanation_stability import importance_stability
def rows(): return [{"feature":"a","importance":2},{"feature":"b","importance":1}]
def test_identical_importance_stable(): assert importance_stability(rows(),rows(),2)["spearman_correlation"]==1
def test_top_overlap_bounded(): assert 0<=importance_stability(rows(),rows(),2)["top_feature_overlap"]<=1
