"""Multiple applicability-domain diagnostics."""
from __future__ import annotations
import numpy as np
from sklearn.covariance import LedoitWolf
from sklearn.neighbors import NearestNeighbors
def mahalanobis(train,test):
 model=LedoitWolf().fit(train); delta=np.asarray(test)-model.location_; return np.sqrt(np.einsum("ij,jk,ik->i",delta,model.precision_,delta))
def knn_distance(train,test,k=5): return NearestNeighbors(n_neighbors=min(k,len(train))).fit(train).kneighbors(test)[0].mean(1)
def local_density(train,test,k=5): return 1/(knn_distance(train,test,k)+1e-12)
def leverage(train,test):
 x=np.c_[np.ones(len(train)),train]; xt=np.c_[np.ones(len(test)),test]; inv=np.linalg.pinv(x.T@x); return np.einsum("ij,jk,ik->i",xt,inv,xt)
