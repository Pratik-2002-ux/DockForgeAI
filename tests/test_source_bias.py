from src.validation.core import SOURCE_FIELDS, load_features
def test_source_fields_are_not_molecular_features():
    _,_,index=load_features("descriptor"); assert not set(SOURCE_FIELDS)&set(index.columns)
def test_source_fields_do_not_include_label(): assert "binary_label" not in SOURCE_FIELDS
