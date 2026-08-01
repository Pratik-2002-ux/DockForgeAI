from __future__ import annotations
import numpy as np
from sklearn.metrics import pairwise_distances
def nearest_neighbors(train,test,train_index,k=5):
 similarities=1-pairwise_distances(test,train,metric="jaccard"); order=np.argsort(-similarities,axis=1)[:,:k]
 return [[{"neighbor_compound_id":train_index.iloc[j].compound_id,"neighbor_label":train_index.iloc[j].binary_label,"tanimoto_similarity":float(similarities[i,j])} for j in indices] for i,indices in enumerate(order)]
