"""Prospective Track 1 inference orchestration."""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from .applicability_domain import applicability_domain
from .feature_rebuilder import rebuild_morgan_r3, schema_checksum, validate_feature_schema
from .model_loader import file_sha256, load_frozen_model
from .ranking import rank_predictions, review_pools
from .uncertainty import uncertainty_columns


REQUIRED_PREPARED = {
    "compound_id", "canonical_isomeric_smiles", "supplier_id", "catalogue_id", "availability",
    "inchikey", "inchikey_connectivity", "murcko_scaffold_smiles", "molecular_weight",
    "exact_molecular_weight", "logp", "tpsa", "hbd", "hba", "rotatable_bonds", "ring_count",
    "aromatic_ring_count", "fraction_csp3", "formal_charge", "qed", "pains_flag", "brenk_flag",
    "extreme_molecular_weight_flag", "extreme_logp_flag", "extreme_tpsa_flag",
    "excessive_rotatable_bonds_flag", "unusual_formal_charge_flag", "structural_review_flags",
}


def _atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    with tempfile.NamedTemporaryFile("w", suffix=".csv", dir=path.parent, delete=False, newline="", encoding="utf-8") as handle:
        temporary = Path(handle.name)
        frame.to_csv(handle, index=False)
    os.replace(temporary, path)


def _atomic_json(value: object, path: Path) -> None:
    with tempfile.NamedTemporaryFile("w", suffix=".json", dir=path.parent, delete=False, encoding="utf-8") as handle:
        temporary = Path(handle.name)
        json.dump(value, handle, indent=2)
        handle.write("\n")
    os.replace(temporary, path)


def _atomic_text(value: str, path: Path) -> None:
    with tempfile.NamedTemporaryFile("w", suffix=".md", dir=path.parent, delete=False, encoding="utf-8") as handle:
        temporary = Path(handle.name)
        handle.write(value)
    os.replace(temporary, path)


def _display_path(path: Path, root: Path) -> str:
    return str(path.relative_to(root)) if path.is_relative_to(root) else str(path)


def _validate_prepared(input_path: Path) -> tuple[pd.DataFrame, dict, Path]:
    manifest_path = input_path.with_name("preparation_manifest.json")
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Preparation manifest unavailable: {manifest_path}")
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("screening_ready") is not True:
        raise ValueError("Prepared library is not screening_ready")
    expected = manifest.get("output_checksums", {}).get(input_path.name)
    if not expected or file_sha256(input_path) != expected:
        raise ValueError("Prepared-library checksum mismatch")
    frame = pd.read_csv(input_path)
    missing = sorted(REQUIRED_PREPARED - set(frame.columns))
    if missing:
        raise ValueError(f"Prepared library missing required columns: {missing}")
    if len(frame) != manifest.get("final_usable_count"):
        raise ValueError("Prepared-library row count disagrees with preparation manifest")
    required_nonempty = ["compound_id", "canonical_isomeric_smiles", "supplier_id", "catalogue_id", "availability", "inchikey", "inchikey_connectivity"]
    if frame[required_nonempty].isna().any().any() or frame[required_nonempty].astype(str).eq("").any().any():
        raise ValueError("Prepared library has missing structure or supplier provenance")
    return frame, manifest, manifest_path


def run_screening(input_path: Path, model_path: Path, output_dir: Path, root: Path, overwrite: bool = False) -> dict:
    output_names = [
        "prospective_predictions.csv", "in_domain_candidates.csv", "borderline_candidates.csv",
        "out_of_domain_candidates.csv", "high_uncertainty_candidates.csv", "structural_alert_review.csv",
        "top_100_review_pool.csv", "top_30_diverse_review_pool.csv", "feature_schema_validation.json",
        "prediction_summary.csv", "screening_manifest.json", "screening_schema.json",
        "screening_checksums.json", "screening_summary.md",
    ]
    existing = [output_dir / name for name in output_names if (output_dir / name).exists()]
    if existing and not overwrite:
        raise FileExistsError(f"Outputs already exist; use --overwrite: {existing[0]}")
    output_dir.mkdir(parents=True, exist_ok=True)
    prepared, preparation, preparation_manifest_path = _validate_prepared(input_path)
    model, metadata, metadata_paths = load_frozen_model(model_path, root)
    training_manifest = metadata["training_manifest"]
    dataset_path = root / "data/model_v2/regression/dataset/ki_regression_dataset.csv"
    if file_sha256(dataset_path) != training_manifest.get("dataset_sha256"):
        raise ValueError("Training dataset checksum mismatch")
    training = pd.read_csv(dataset_path)
    fitted = training[training["dataset_split"].isin(["train", "validation"])].reset_index(drop=True)
    feature_artifact = root / "data/model_v2/regression/features/morgan_r3.npz"
    persisted = np.load(feature_artifact)["X"]
    rebuilt_training_all = rebuild_morgan_r3(training["canonical_isomeric_smiles"].tolist())
    if persisted.shape != rebuilt_training_all.shape or not np.array_equal(persisted, rebuilt_training_all):
        raise ValueError("Rebuilt feature order disagrees with the frozen Track 1 feature artifact")
    schema_report = validate_feature_schema(model, training_manifest.get("feature_set"))
    schema_report.update({"training_feature_artifact": str(feature_artifact.relative_to(root)), "training_feature_artifact_checksum": file_sha256(feature_artifact), "training_feature_rows_validated": len(training)})
    if not schema_report["matches_frozen_model"]:
        raise ValueError("Prospective feature schema does not match frozen model")
    _atomic_json(schema_report, output_dir / "feature_schema_validation.json")
    print(f"Rebuilding frozen features for {len(prepared):,} compounds", flush=True)
    features = rebuild_morgan_r3(prepared["canonical_isomeric_smiles"].tolist())
    predictions = model.predict(features)
    threshold = float(metadata["best_model"]["applicability_domain"]["threshold_max_tanimoto"])
    training_features = rebuilt_training_all[training["dataset_split"].isin(["train", "validation"]).to_numpy()]
    domain = applicability_domain(features, training_features, fitted["compound_id"], fitted["median_pki"], threshold)
    uncertainty, uncertainty_meta = uncertainty_columns(predictions, metadata)
    result = prepared.copy()
    result["predicted_pki"] = predictions
    result["predicted_ki_nm"] = 10 ** (9 - predictions)
    result["model_id"] = training_manifest["model_id"]
    result["model_family"] = training_manifest["model_family"]
    result["feature_representation"] = training_manifest["feature_set"]
    for name, values in uncertainty.items():
        result[name] = values
    for name, values in domain.items():
        result[name] = values
    result = rank_predictions(result)
    top100, top30 = review_pools(result)
    subsets = {
        "prospective_predictions.csv": result,
        "in_domain_candidates.csv": result[result.applicability_domain_status.eq("in_domain")],
        "borderline_candidates.csv": result[result.applicability_domain_status.eq("borderline")],
        "out_of_domain_candidates.csv": result[result.applicability_domain_status.eq("out_of_domain")],
        "high_uncertainty_candidates.csv": result.iloc[0:0],
        "structural_alert_review.csv": result[result.pains_flag.astype(bool) | result.brenk_flag.astype(bool) | result.structural_review_flags.fillna("").ne("")],
        "top_100_review_pool.csv": top100,
        "top_30_diverse_review_pool.csv": top30,
    }
    for name, frame in subsets.items():
        _atomic_csv(frame, output_dir / name)
    counts = result.applicability_domain_status.value_counts()
    summary = pd.DataFrame([
        {"metric": "prediction_count", "value": len(result)},
        {"metric": "in_domain_count", "value": int(counts.get("in_domain", 0))},
        {"metric": "borderline_count", "value": int(counts.get("borderline", 0))},
        {"metric": "out_of_domain_count", "value": int(counts.get("out_of_domain", 0))},
        {"metric": "predicted_pki_min", "value": float(result.predicted_pki.min())},
        {"metric": "predicted_pki_max", "value": float(result.predicted_pki.max())},
        {"metric": "predicted_ki_nm_min", "value": float(result.predicted_ki_nm.min())},
        {"metric": "predicted_ki_nm_max", "value": float(result.predicted_ki_nm.max())},
    ])
    _atomic_csv(summary, output_dir / "prediction_summary.csv")
    data_outputs = [output_dir / name for name in [*subsets, "feature_schema_validation.json", "prediction_summary.csv"]]
    output_checksums = {path.name: file_sha256(path) for path in data_outputs}
    docking_ready = bool(len(result) == len(prepared) and schema_report["matches_frozen_model"] and int(counts.get("in_domain", 0)) > 0 and len(top30) == min(30, len(result)) and all(file_sha256(path) == output_checksums[path.name] for path in data_outputs))
    warnings = [
        "Predictions are model estimates for research review and are not inhibitor or activity claims.",
        "Track 1 records only a binary applicability-domain rule; no borderline band was invented.",
        "The deployable artifact is one estimator, so per-compound ensemble standard deviation is unavailable.",
        "Conformal intervals use only the half-width persisted in frozen Track 1 metadata.",
        "Availability is from the supplier snapshot and does not guarantee current stock.",
        "No docking was run and no experimental shortlist was frozen.",
    ]
    manifest = {
        "schema_version": 1, "study_id": preparation["study_id"],
        "prepared_library_path": _display_path(input_path, root), "prepared_library_checksum": file_sha256(input_path),
        "preparation_manifest_path": _display_path(preparation_manifest_path, root),
        "frozen_model_path": _display_path(model_path, root), "frozen_model_checksum": file_sha256(model_path),
        "model_card_path": str(metadata_paths["model_card"].relative_to(root)),
        "training_manifest_path": str(metadata_paths["training_manifest"].relative_to(root)),
        "best_model_metadata_path": str(metadata_paths["best_model"].relative_to(root)),
        "training_dataset_checksum": training_manifest["dataset_sha256"], "feature_schema_checksum": schema_checksum(),
        "training_code_git_commit": training_manifest["git_commit"],
        "current_inference_git_commit": subprocess.run(["git", "rev-parse", "HEAD"], cwd=root, text=True, capture_output=True, check=True).stdout.strip(),
        "input_count": len(prepared), "prediction_count": len(result),
        "in_domain_count": int(counts.get("in_domain", 0)), "borderline_count": int(counts.get("borderline", 0)), "out_of_domain_count": int(counts.get("out_of_domain", 0)),
        "uncertainty_available": uncertainty_meta["available"], "uncertainty_method": uncertainty_meta["method"], "high_uncertainty_count": 0,
        "pains_flag_count": int(result.pains_flag.astype(bool).sum()), "brenk_flag_count": int(result.brenk_flag.astype(bool).sum()),
        "top_100_count": len(top100), "top_30_diverse_count": len(top30),
        "output_checksums": output_checksums, "generated_timestamp": datetime.now(timezone.utc).isoformat(),
        "warnings": warnings, "docking_ready": docking_ready,
    }
    schema = {"prediction_columns": list(result.columns), "output_columns": {name: list(frame.columns) for name, frame in subsets.items()}}
    _atomic_json(schema, output_dir / "screening_schema.json")
    _atomic_json(manifest, output_dir / "screening_manifest.json")
    markdown = f"# Prospective Track 1 ML screening\n\n- Predictions: {len(result):,}\n- In domain: {manifest['in_domain_count']:,}\n- Out of domain: {manifest['out_of_domain_count']:,}\n- Uncertainty: {manifest['uncertainty_method']}\n- Docking ready: {docking_ready}\n\nNo docking was run and no experimental shortlist was frozen.\n"
    _atomic_text(markdown, output_dir / "screening_summary.md")
    checksum_paths = [*data_outputs, output_dir / "screening_schema.json", output_dir / "screening_manifest.json", output_dir / "screening_summary.md"]
    _atomic_json({path.name: file_sha256(path) for path in checksum_paths}, output_dir / "screening_checksums.json")
    return manifest
