import json,pickle,pytest
from src.models_v2.export import export_pickle,export_metadata,export_onnx,export_torchscript
from src.models_v2.reporting import write_publication_bundle
def test_pickle_export(tmp_path): assert pickle.loads(export_pickle({"a":1},tmp_path/"x.pkl").read_bytes())=={"a":1}
def test_metadata_export(tmp_path): assert json.loads(export_metadata({"a":1},tmp_path/"x.json").read_text())["a"]==1
@pytest.mark.parametrize("fn",[export_onnx,export_torchscript])
def test_optional_export_fails_clearly(fn):
 with pytest.raises(RuntimeError): fn()
def test_publication_bundle(tmp_path): assert len(write_publication_bundle(tmp_path,{}))==5
