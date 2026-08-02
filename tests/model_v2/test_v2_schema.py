import pandas as pd,pytest
from src.models_v2.schema import validate_dataset,V2SchemaError
def valid(): return pd.DataFrame({"compound_id":["a"],"canonical_isomeric_smiles":["CCO"],"binary_label":["active"],"source_database":["ChEMBL"],"source_assay_id":["x"],"measurement_type":["IC50"],"standard_value_nm":[1.],"assay_date":["2020-01-01"],"murcko_scaffold_smiles":["CCO"],"record_weight":[1.]})
def test_valid_schema(): assert validate_dataset(valid())
@pytest.mark.parametrize("label",["decoy","unknown","positive","negative"])
def test_nonexperimental_labels_rejected(label):
 d=valid(); d.binary_label=label
 with pytest.raises(V2SchemaError): validate_dataset(d)
@pytest.mark.parametrize("column",valid().columns)
def test_required_columns_enforced(column):
 with pytest.raises(V2SchemaError): validate_dataset(valid().drop(columns=column))
