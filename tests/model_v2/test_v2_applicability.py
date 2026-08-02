import numpy as np,pytest
from src.models_v2.applicability import mahalanobis,knn_distance,local_density,leverage
TRAIN=np.array([[0.,0.],[1,0],[0,1],[1,1]])
@pytest.mark.parametrize("fn",[mahalanobis,knn_distance,local_density,leverage])
def test_domain_metric_shape(fn): assert fn(TRAIN,np.array([[.5,.5]])).shape==(1,)
def test_knn_distance_nonnegative(): assert (knn_distance(TRAIN,TRAIN)>=0).all()
def test_density_positive(): assert (local_density(TRAIN,TRAIN)>0).all()
