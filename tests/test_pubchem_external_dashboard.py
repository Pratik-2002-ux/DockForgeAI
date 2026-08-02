from app.components.pubchem_external import DOWNLOAD_KEYS,build_bundle
def test_widget_keys_are_unique(): assert len(DOWNLOAD_KEYS)==len(set(DOWNLOAD_KEYS.values()))
def test_empty_bundle_is_zip(tmp_path):
 (tmp_path/"data/processed/pubchem_external_validation").mkdir(parents=True); assert build_bundle(tmp_path).startswith(b"PK")
