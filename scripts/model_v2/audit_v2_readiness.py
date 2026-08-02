#!/usr/bin/env python3
"""Audit approved local v2 inputs without training or network access."""
from __future__ import annotations

import argparse
import importlib.metadata
import importlib.util
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml
from rdkit import RDLogger

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.datasets.activity_curator import curate_activities  # noqa: E402
from src.datasets.chembl_client import ChEMBLClient  # noqa: E402


PACKAGES = (
    "numpy", "pandas", "scipy", "scikit-learn", "joblib", "pyarrow",
    "matplotlib", "streamlit", "altair", "shap", "rdkit", "pyyaml",
    "tqdm", "pydantic", "reportlab", "pypdf", "xgboost", "lightgbm",
    "catboost", "torch", "torch-geometric", "onnx", "onnxruntime",
    "skl2onnx",
)


def _available(distribution: str) -> dict[str, object]:
    module = {
        "scikit-learn": "sklearn", "pyyaml": "yaml",
        "torch-geometric": "torch_geometric",
    }.get(distribution, distribution.replace("-", "_"))
    present = importlib.util.find_spec(module) is not None
    try:
        version = importlib.metadata.version(distribution) if present else None
    except importlib.metadata.PackageNotFoundError:
        version = None
    return {"available": present, "version": version}


def audit() -> dict[str, object]:
    settings = yaml.safe_load((ROOT / "config/settings.yaml").read_text())
    config = settings["xai_dataset"]
    chembl_dir = ROOT / settings["directories"]["xai_chembl_raw"]
    pages = sorted(chembl_dir.glob("activity_page_*.json"))
    RDLogger.DisableLog("rdApp.*")
    client = ChEMBLClient(
        chembl_dir,
        config["request_timeout_seconds"],
        config["request_retries"],
        config["request_page_size"],
    )
    raw = client.fetch_activities(config["target_chembl_id"], False, True)
    curated = curate_activities(raw, config)
    locked_path = ROOT / "data/processed/pubchem_external_validation/dataset/pubchem_external_validation.csv"
    bindingdb_exports = sorted(ROOT.glob("data/**/*bindingdb*.csv"))
    eligible_sources = ["ChEMBL"] if pages else []
    preliminary = {
        "raw_records": len(raw),
        **{key: len(value) for key, value in curated.items()},
    }
    blockers = []
    if not bindingdb_exports:
        blockers.append("No approved local BindingDB training export is present.")
    blockers.append("The local PubChem data are the locked external benchmark and are excluded from training.")
    if preliminary["inactive"] < 200:
        blockers.append("Fewer than 200 curated experimental inactives are available before locked-overlap exclusion.")
    else:
        blockers.append("Only 227 preliminary ChEMBL inactives exist; locked-overlap curation may reduce the final class below 200.")
    blockers.append("A scientifically suitable multi-source training dataset cannot be constructed from current approved local files.")
    return {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "offline": True,
        "training_allowed": False,
        "status": "readiness_only",
        "approved_training_sources": eligible_sources,
        "locked_external": {
            "path": str(locked_path.relative_to(ROOT)),
            "exists": locked_path.is_file(),
            "training_use_forbidden": True,
        },
        "bindingdb_training_exports": [str(p.relative_to(ROOT)) for p in bindingdb_exports],
        "chembl_preliminary_audit": preliminary,
        "dependencies": {name: _available(name) for name in PACKAGES},
        "blockers": blockers,
        "required_next_input": "An explicitly approved, local experimental inactive dataset with complete assay provenance from an independent source.",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="data/model_v2/readiness/readiness_report.json")
    args = parser.parse_args()
    report = audit()
    output = ROOT / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

