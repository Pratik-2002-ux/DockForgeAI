import pandas as pd,pytest
from src.models_v2.splits import scaffold_split,assert_no_scaffold_leakage,time_split,leave_group_out
def frame(): return pd.DataFrame({"murcko_scaffold_smiles":[f"s{i//2}" for i in range(30)],"assay_date":pd.date_range("2020-01-01",periods=30),"source_database":["a"]*15+["b"]*15})
@pytest.mark.parametrize("seed",[1,2,3,2026,2030])
def test_scaffold_split_never_leaks(seed):
 d=frame(); split=scaffold_split(d,seed); assert assert_no_scaffold_leakage(d,split)
def test_scaffold_split_deterministic(): assert scaffold_split(frame(),7).equals(scaffold_split(frame(),7))
def test_time_split_has_all_partitions(): assert set(time_split(frame()))=={"train","validation","test"}
def test_leave_source_out_has_two_groups(): assert len(leave_group_out(frame(),"source_database"))==2
