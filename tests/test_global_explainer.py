from src.explainability.global_explainer import model_importance
from sklearn.linear_model import LogisticRegression
def model(): return LogisticRegression().fit([[0,0],[1,1],[0,1],[1,0]],[0,1,0,1])
def test_global_importance_sorted():
 rows=model_importance(model(),["a","b"]);assert rows[0]["importance"]>=rows[1]["importance"]
def test_global_importance_has_direction(): assert "direction" in model_importance(model(),["a","b"])[0]
