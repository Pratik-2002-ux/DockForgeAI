"""Validation and loading for the frozen Track 1 regression model."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import joblib


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_frozen_model(model_path: Path, root: Path):
    paths = {
        "model_card": root / "data/model_v2/regression/models/model_card.json",
        "training_manifest": root / "data/model_v2/regression/models/training_manifest.json",
        "best_model": root / "data/model_v2/regression/evaluation/best_model.json",
    }
    missing = [str(path) for path in [model_path, *paths.values()] if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing frozen model artifact(s): " + ", ".join(missing))
    metadata = {name: json.loads(path.read_text()) for name, path in paths.items()}
    training = metadata["training_manifest"]
    card_metrics = metadata["model_card"].get("metrics", {})
    best = metadata["best_model"]
    if file_sha256(model_path) != training.get("model_sha256"):
        raise ValueError("Frozen model checksum mismatch")
    comparisons = (
        (training.get("model_family"), best.get("best_model"), "model family"),
        (training.get("feature_set"), best.get("best_feature_set"), "feature representation"),
        (training.get("model_family"), card_metrics.get("best_model"), "model-card family"),
        (training.get("feature_set"), card_metrics.get("best_feature_set"), "model-card features"),
    )
    for left, right, label in comparisons:
        if left != right:
            raise ValueError(f"Frozen metadata disagree on {label}: {left!r} != {right!r}")
    model = joblib.load(model_path)
    if getattr(model, "n_features_in_", None) != 2048 or training.get("feature_set") != "morgan_r3":
        raise ValueError("Frozen estimator feature order/size does not match Morgan radius 3, 2048 bits")
    return model, metadata, paths
