import json
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[2]

def test_binary_training_is_disabled():
 config=json.loads((ROOT/"data/model_v2/binary/future_training_config.json").read_text());assert config["training_allowed"] is False and config["enabled"] is False
def test_readiness_gate_remains_blocked():
 report=json.loads((ROOT/"data/model_v2/readiness/readiness_report.json").read_text());assert report["training_allowed"] is False;assert report["source_shortcut_risk"] is True
def test_regression_uses_only_ki():
 manifest=json.loads((ROOT/"data/model_v2/regression/dataset/ki_dataset_manifest.json").read_text());assert manifest["endpoint"]=="Ki" and manifest["ordinary_regression_relations"]==["="]
def test_locked_external_not_used_for_selection():
 best=json.loads((ROOT/"data/model_v2/regression/evaluation/best_model.json").read_text());assert best["selection_uses_external"] is False
def test_track_two_not_in_prospective_module():
 text="".join(p.read_text() for p in (ROOT/"src/prospective_screening").glob("*.py"));assert "binary_label" not in text and "classification" not in text
@pytest.mark.parametrize("artifact",["ki_dataset_manifest.json","ki_dataset_schema.json","ki_dataset_checksums.json"])
def test_regression_dataset_manifests(artifact):assert (ROOT/"data/model_v2/regression/dataset"/artifact).is_file()
@pytest.mark.parametrize("artifact",["training_manifest.json","model_card.json"])
def test_model_manifests(artifact):assert (ROOT/"data/model_v2/regression/models"/artifact).is_file()
