from src.models_v2.data.split_builder import _group_split,leakage_count
def assign_scaffold_split(frame):
 out=frame.copy(); out["dataset_split"]=_group_split(out.murcko_scaffold_smiles); return out
