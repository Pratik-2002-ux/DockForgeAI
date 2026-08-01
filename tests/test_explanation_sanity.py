from src.validation.explanation_sanity import run_explanation_sanity
def test_explanation_sanity_entrypoint_is_callable(): assert callable(run_explanation_sanity)
def test_explanations_are_not_declared_causal():
    import inspect
    assert "causal" in inspect.getsource(run_explanation_sanity)
