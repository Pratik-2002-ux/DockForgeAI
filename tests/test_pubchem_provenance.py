from src.external_validation.pubchem_provenance import checksum,write_manifest
def test_checksum_is_deterministic(tmp_path):
 p=tmp_path/"x"; p.write_text("x"); assert checksum(p)==checksum(p)
def test_manifest_marks_model_locked(tmp_path):
 p=tmp_path/"x"; p.write_text("x"); assert write_manifest(tmp_path/"m.json",[p])["locked_model_retrained"] is False
