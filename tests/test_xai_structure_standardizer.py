import pytest
from src.datasets.structure_standardizer import StructureStandardizationError, standardize_structure

def test_standardization_removes_counterion_and_describes() -> None:
    result=standardize_structure("CCO.[Na+]")
    assert result["canonical_isomeric_smiles"]=="CCO" and result["heavy_atom_count"]==3 and result["inchikey"]
def test_invalid_structure_rejected() -> None:
    with pytest.raises(StructureStandardizationError): standardize_structure("not[smiles")
