from pathlib import Path
from app.components.xai_validation import external_template,validation_bundle
def test_external_template_columns(): assert external_template().startswith(b"compound_id,canonical_smiles,true_label")
def test_validation_bundle_is_zip(tmp_path):
    (tmp_path/"data/processed/xai_validation").mkdir(parents=True); (tmp_path/"data/processed/xai_validation/a.csv").write_text("a\n1\n"); assert validation_bundle(tmp_path).startswith(b"PK")
