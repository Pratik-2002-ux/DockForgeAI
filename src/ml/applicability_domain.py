from __future__ import annotations
import numpy as np
from sklearn.metrics import pairwise_distances
from sklearn.preprocessing import StandardScaler
def applicability_domain(train_fp,test_fp,train_desc,test_desc,train_ids,quantile=.05):
 sim=1-pairwise_distances(test_fp,train_fp,metric="jaccard"); order=np.argsort(-sim,axis=1)[:,:5]; scaler=StandardScaler().fit(train_desc); td=scaler.transform(train_desc); qd=scaler.transform(test_desc); distances=pairwise_distances(qd,td).min(axis=1); train_dist=pairwise_distances(td,td); np.fill_diagonal(train_dist,np.nan); train_sim=1-pairwise_distances(train_fp,train_fp,metric="jaccard"); np.fill_diagonal(train_sim,np.nan); threshold_sim=float(np.quantile(np.nanmax(train_sim,axis=1),quantile)); threshold_dist=float(np.quantile(np.nanmin(train_dist,axis=1),.95))
 return [{"nearest_training_compound_id":train_ids[idx[0]],"maximum_tanimoto_similarity":float(sim[i,idx[0]]),"mean_top5_similarity":float(sim[i,idx].mean()),"descriptor_distance":float(distances[i]),"in_domain":bool(sim[i,idx[0]]>=threshold_sim and distances[i]<=threshold_dist),"domain_confidence_score":float(max(0,min(1,sim[i,idx[0]]))),"warning_reason":"" if sim[i,idx[0]]>=threshold_sim and distances[i]<=threshold_dist else "low similarity or high descriptor distance"} for i,idx in enumerate(order)]
