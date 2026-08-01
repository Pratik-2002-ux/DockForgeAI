import pandas as pd, pytest
from src.validation.external_dataset_validator import validate_external_dataset
def test_external_schema_validation(tmp_path):
    p=tmp_path/"input.csv"; pd.DataFrame({"id":[1]}).to_csv(p,index=False)
    with pytest.raises(ValueError): validate_external_dataset(p,"id","smiles",tmp_path/"out")
def test_external_validation_source_has_no_fit_call():
    import inspect
    source=inspect.getsource(validate_external_dataset); assert ".fit(" not in source and "threshold_optimized\":False" in source
