from src.validation.calibration_robustness import maximum_calibration_error,calibration_slope_intercept
def test_perfect_calibration_has_small_error(): assert maximum_calibration_error([0,1],[0,1])<1e-8
def test_slope_intercept_are_finite():
    import numpy as np
    s,i=calibration_slope_intercept([0,0,1,1],[.1,.2,.8,.9]); assert np.isfinite([s,i]).all()
