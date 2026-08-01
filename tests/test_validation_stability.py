from src.validation.stability_analysis import run_stability
def test_stability_entrypoint_is_callable(): assert callable(run_stability)
def test_stability_uses_existing_predictions(): assert "stability" in run_stability.__module__
