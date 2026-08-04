from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.prospective_screening.applicability_domain import applicability_domain
from src.prospective_screening.feature_rebuilder import rebuild_morgan_r3, validate_feature_schema
from src.prospective_screening.inference import run_screening
from src.prospective_screening.model_loader import file_sha256, load_frozen_model
from src.prospective_screening.ranking import rank_predictions, review_pools
from src.prospective_screening.uncertainty import uncertainty_columns


ROOT = Path(__file__).resolve().parents[1]
MODEL = ROOT / "data/model_v2/regression/models/best_ki_regressor.joblib"
SOURCE = ROOT / "data/processed/prospective_studies/HIVPR_2026_001/prepared_library.csv"


def prepared_fixture(tmp_path: Path, rows: int = 4) -> Path:
    frame = pd.read_csv(SOURCE).head(rows)
    training = pd.read_csv(ROOT / "data/model_v2/regression/dataset/ki_regression_dataset.csv")
    frame.loc[frame.index[0], "canonical_isomeric_smiles"] = training.loc[0, "canonical_isomeric_smiles"]
    path = tmp_path / "prepared_library.csv"
    frame.to_csv(path, index=False)
    checksum = hashlib.sha256(path.read_bytes()).hexdigest()
    (tmp_path / "preparation_manifest.json").write_text(json.dumps({"study_id": "TEST", "screening_ready": True, "final_usable_count": rows, "output_checksums": {"prepared_library.csv": checksum}}))
    return path


@pytest.fixture()
def screened(tmp_path):
    source = prepared_fixture(tmp_path)
    output = tmp_path / "screening"
    manifest = run_screening(source, MODEL, output, ROOT)
    return output, manifest


def test_frozen_model_loading():
    model, metadata, _ = load_frozen_model(MODEL, ROOT)
    assert model.n_features_in_ == 2048
    assert metadata["training_manifest"]["model_family"] == "random_forest"


def test_model_checksum_validation(tmp_path):
    broken = tmp_path / "model.joblib"
    broken.write_bytes(MODEL.read_bytes() + b"broken")
    with pytest.raises(ValueError, match="checksum"):
        load_frozen_model(broken, ROOT)


def test_manifest_consistency():
    _, metadata, _ = load_frozen_model(MODEL, ROOT)
    assert metadata["training_manifest"]["feature_set"] == metadata["best_model"]["best_feature_set"]


def test_missing_model_artifact(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_frozen_model(tmp_path / "missing.joblib", ROOT)


def test_feature_order_mismatch():
    class Wrong:
        n_features_in_ = 9
    assert validate_feature_schema(Wrong(), "morgan_r3")["matches_frozen_model"] is False


def test_prepared_library_checksum_mismatch(tmp_path):
    source = prepared_fixture(tmp_path)
    source.write_text(source.read_text() + "\n")
    with pytest.raises(ValueError, match="checksum"):
        run_screening(source, MODEL, tmp_path / "out", ROOT)


def test_screening_ready_false(tmp_path):
    source = prepared_fixture(tmp_path)
    manifest_path = tmp_path / "preparation_manifest.json"
    manifest = json.loads(manifest_path.read_text()); manifest["screening_ready"] = False
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="screening_ready"):
        run_screening(source, MODEL, tmp_path / "out", ROOT)


def test_exact_prediction_count(screened):
    output, manifest = screened
    assert len(pd.read_csv(output / "prospective_predictions.csv")) == manifest["input_count"] == 4


def test_pki_to_ki_conversion(screened):
    predictions = pd.read_csv(screened[0] / "prospective_predictions.csv")
    assert np.allclose(predictions.predicted_ki_nm, 10 ** (9 - predictions.predicted_pki))


def test_applicability_domain_threshold():
    fp = rebuild_morgan_r3(["CCO"])
    domain = applicability_domain(fp, fp, ["x"], [7.0], 0.3)
    assert domain["applicability_domain_status"].tolist() == ["in_domain"]
    assert domain["maximum_training_tanimoto"][0] == 1.0


def test_uncertainty_unavailable_handling():
    values, meta = uncertainty_columns(np.array([7.0]), {"best_model": {"uncertainty": {}}})
    assert meta["available"] is False
    assert values["uncertainty_status"][0] == "unavailable"


def test_ranking_determinism(screened):
    frame = pd.read_csv(screened[0] / "prospective_predictions.csv")
    reranked = rank_predictions(frame.drop(columns=["predicted_pki_rank", "uncertainty_rank", "domain_rank", "diversity_rank", "review_priority_score", "review_priority_rank"]))
    assert reranked.review_priority_rank.tolist() == frame.review_priority_rank.tolist()


def test_structural_alert_retention(screened):
    predictions = pd.read_csv(screened[0] / "prospective_predictions.csv")
    assert {"pains_flag", "brenk_flag", "structural_review_flags"} <= set(predictions)


def test_top_100_generation(screened):
    assert len(pd.read_csv(screened[0] / "top_100_review_pool.csv")) == 4


def test_diverse_top_30_generation(screened):
    predictions = pd.read_csv(screened[0] / "prospective_predictions.csv")
    _, diverse = review_pools(predictions)
    assert len(diverse) == 4


def test_overwrite_behavior(screened):
    output, _ = screened
    with pytest.raises(FileExistsError):
        run_screening(output.parent / "prepared_library.csv", MODEL, output, ROOT)


def test_screening_manifest_generation(screened):
    output, manifest = screened
    assert json.loads((output / "screening_manifest.json").read_text())["prediction_count"] == manifest["prediction_count"]


def test_checksum_generation(screened):
    output, _ = screened
    checks = json.loads((output / "screening_checksums.json").read_text())
    assert all(file_sha256(output / name) == checksum for name, checksum in checks.items())


def test_docking_ready_logic(screened):
    assert screened[1]["docking_ready"] is True
    assert screened[1]["top_30_diverse_count"] == 4
