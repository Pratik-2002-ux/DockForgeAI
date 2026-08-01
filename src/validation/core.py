"""Shared, deterministic validation helpers."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.metrics import (
    accuracy_score, average_precision_score, balanced_accuracy_score,
    brier_score_loss, confusion_matrix, log_loss, matthews_corrcoef,
    roc_auc_score,
)

ROOT = Path(__file__).resolve().parents[2]
DATASET = ROOT / "data/processed/xai_dataset/hivpr_xai_2000_compounds.csv"
FEATURES = ROOT / "data/processed/xai_features"
PREDICTIONS = ROOT / "data/processed/xai_predictions"
BENCHMARKS = ROOT / "data/processed/xai_benchmarks"
OUTPUT = ROOT / "data/processed/xai_validation"
BEST_MODEL_ID = "logistic_regression__fingerprint__seed2026"
SOURCE_FIELDS = [
    "source_database", "label_source", "source_target_id", "source_assay_id",
    "source_document_id", "provenance_reference", "retrieval_date", "source_release",
]

def ensure_dirs(root: Path = OUTPUT) -> None:
    for name in ("source_bias", "permutation", "baselines", "similarity", "scaffolds",
                 "matched_subsets", "reweighting", "thresholds", "decoy_sensitivity",
                 "calibration", "stability", "explanation_sanity", "probability_shift",
                 "figures", "reports"):
        (root / name).mkdir(parents=True, exist_ok=True)
    (ROOT / "logs/xai_validation").mkdir(parents=True, exist_ok=True)

def load_data() -> pd.DataFrame:
    return pd.read_csv(DATASET)

def load_features(name: str = "fingerprint"):
    index = pd.read_csv(FEATURES / "feature_index.csv")
    if name == "descriptor":
        x = pd.read_parquet(FEATURES / "descriptor_features.parquet").to_numpy(float)
    elif name == "fingerprint":
        x = sparse.load_npz(FEATURES / "morgan_features.npz")
    elif name == "fused":
        x = sparse.load_npz(FEATURES / "fused_features.npz")
    else:
        raise ValueError(f"Unknown feature family: {name}")
    return x, (index.binary_label == "active").astype(int).to_numpy(), index

def best_predictions(split: str = "test") -> pd.DataFrame:
    frame = pd.read_csv(PREDICTIONS / f"{split}_predictions.csv")
    return frame.loc[frame.model_id.eq(BEST_MODEL_ID)].reset_index(drop=True)

def safe_metrics(y, p) -> dict:
    y = np.asarray(y, int); p = np.clip(np.asarray(p, float), 1e-12, 1 - 1e-12)
    pred = (p >= .5).astype(int)
    result = {
        "sample_count": int(len(y)), "positive_count": int(y.sum()),
        "accuracy": float(accuracy_score(y, pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y, pred)) if len(np.unique(y)) == 2 else float(accuracy_score(y, pred)),
        "mcc": float(matthews_corrcoef(y, pred)) if len(np.unique(y)) == 2 else 0.0,
        "brier_score": float(brier_score_loss(y, p)),
        "log_loss": float(log_loss(y, p, labels=[0, 1])),
        "mean_confidence": float(np.maximum(p, 1-p).mean()),
        "error_rate": float(np.mean(pred != y)),
    }
    result["roc_auc"] = float(roc_auc_score(y, p)) if len(np.unique(y)) == 2 else None
    result["pr_auc"] = float(average_precision_score(y, p)) if y.sum() else None
    tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()
    result.update(tn=int(tn), fp=int(fp), fn=int(fn), tp=int(tp))
    return result

def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n")

def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def jaccard_similarity(query, reference):
    query = sparse.csr_matrix(query).astype(bool).astype(np.uint8)
    reference = sparse.csr_matrix(reference).astype(bool).astype(np.uint8)
    intersection = (query @ reference.T).toarray().astype(float)
    union = np.asarray(query.sum(1)) + np.asarray(reference.sum(1)).T - intersection
    return np.divide(intersection, union, out=np.zeros_like(intersection), where=union != 0)
