import numpy as np
from scipy import sparse
from src.validation.core import jaccard_similarity
def test_tanimoto_identity(): assert jaccard_similarity(sparse.csr_matrix([[1,0,1]]),sparse.csr_matrix([[1,0,1]])).item()==1
def test_tanimoto_disjoint(): assert jaccard_similarity([[1,0]],[[0,1]]).item()==0
