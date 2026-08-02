import numpy as np,pytest
from src.models_v2.calibration import ece,temperature_scale,calibration_metrics
from src.models_v2.uncertainty import ensemble_uncertainty,conformal_threshold,conformal_sets
@pytest.mark.parametrize("bins",[2,5,10,20])
def test_ece_bounded(bins): assert 0<=ece([0,1],[.1,.9],bins)<=1
def test_temperature_positive(): assert temperature_scale([0,1],[.2,.8])[0]>0
def test_brier_present(): assert "brier_score" in calibration_metrics([0,1],[.2,.8])
def test_ensemble_shapes(): assert ensemble_uncertainty([[.1,.9],[.2,.8]])["mean"].shape==(2,)
def test_conformal_threshold_bounded(): assert 0<=conformal_threshold([0,1],[.1,.9])<=1
def test_conformal_sets_count(): assert len(conformal_sets([.1,.9],.2))==2
