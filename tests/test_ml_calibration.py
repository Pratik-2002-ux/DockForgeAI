import numpy as np
from src.ml.calibration import ProbabilityCalibrator,choose_calibration
def test_sigmoid_calibrator_bounds():
 p=ProbabilityCalibrator("sigmoid").fit(np.array([.1,.2,.8,.9]),np.array([0,0,1,1])).predict(np.array([.5]));assert 0<=p[0]<=1
def test_calibration_selection_uses_validation_arrays(): assert choose_calibration(np.array([0,0,1,1]),np.array([.1,.2,.8,.9]))[0] in {"none","sigmoid","isotonic"}
