from src.validation.scaffold_stratification import scaffold_category
def test_scaffold_frequency_categories(): assert [scaffold_category(x) for x in (1,2,4,11)]==["singleton","low","medium","high"]
def test_scaffold_leakage_remains_zero():
    import pandas as pd
    d=pd.read_csv("data/processed/xai_dataset/hivpr_xai_2000_compounds.csv"); assert (d.groupby("scaffold_id").dataset_split.nunique()>1).sum()==0
