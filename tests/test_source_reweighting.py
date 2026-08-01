import pandas as pd
from src.validation.reweighting import source_weights
def test_source_weights_are_clipped():
    d=pd.DataFrame({"source_database":["a"]*9+["b"],"binary_label":["x"]*10}); w=source_weights(d,"source"); assert w.min()>=.1 and w.max()<=10
def test_joint_weights_are_positive():
    d=pd.DataFrame({"source_database":["a","b"],"binary_label":["x","y"]}); assert (source_weights(d)>0).all()
