#!/usr/bin/env python3
"""Build the curated Phase XAI-1 HIVPR active/decoy benchmark dataset."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import yaml
from rdkit import RDLogger

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0, str(ROOT))

from src.datasets.activity_curator import curate_activities  # noqa: E402
from src.datasets.chembl_client import ChEMBLClient  # noqa: E402
from src.datasets.dataset_report import generate_dataset_reports, write_manifest_and_schema, write_sdf, write_table  # noqa: E402
from src.datasets.dataset_validator import validate_output_files, validate_records  # noqa: E402
from src.datasets.decoy_matcher import select_matched_decoys  # noqa: E402
from src.datasets.dude_client import DUDEClient  # noqa: E402
from src.datasets.duplicate_resolver import resolve_duplicates  # noqa: E402
from src.datasets.provenance import retrieval_timestamp, write_provenance  # noqa: E402
from src.datasets.scaffold_splitter import assign_scaffold_splits  # noqa: E402
from src.datasets.structure_standardizer import StructureStandardizationError, standardize_structure  # noqa: E402


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-count", type=int, default=2000); parser.add_argument("--active-count", type=int, default=1000)
    parser.add_argument("--decoy-count", type=int, default=1000); parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--resume", action="store_true"); parser.add_argument("--offline", action="store_true")
    parser.add_argument("--strict", action="store_true"); parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--output-prefix", default="hivpr_xai_2000")
    return parser.parse_args()


def _row_active(record: dict[str, Any], release: str, retrieved: str) -> dict[str, Any]:
    return {
        "compound_id": f"CHEMBL_{record.get('molecule_chembl_id')}", "standardized_compound_id": record["inchikey"], "source_database": "ChEMBL",
        "source_record_id": record["source_record_id"], "source_molecule_id": record.get("molecule_chembl_id"), "source_target_id": record.get("target_chembl_id"),
        "source_assay_id": record.get("assay_chembl_id"), "source_document_id": record.get("document_chembl_id"), "original_smiles": record["original_smiles"],
        **{key: record.get(key) for key in ("canonical_isomeric_smiles", "inchi", "inchikey", "inchikey_connectivity", "molecular_formula", "molecular_weight", "exact_molecular_weight", "logp", "tpsa", "hbd", "hba", "rotatable_bonds", "ring_count", "aromatic_ring_count", "formal_charge", "heavy_atom_count", "fraction_csp3", "qed", "pains_alert_count", "brenk_alert_count", "murcko_scaffold_smiles", "scaffold_id")},
        "activity_type": record.get("standard_type"), "original_activity_value": record.get("standard_value"), "original_activity_units": record.get("standard_units"),
        "standard_activity_value_nm": record.get("standard_activity_value_nm"), "pchembl_value": record.get("activity_median_pchembl"), "activity_relation": record.get("standard_relation"),
        "binary_label": "active", "label_source": "chembl_activity", "activity_record_count": record.get("activity_record_count"),
        "activity_min_pchembl": record.get("activity_min_pchembl"), "activity_max_pchembl": record.get("activity_max_pchembl"), "activity_median_pchembl": record.get("activity_median_pchembl"),
        "activity_conflict_flag": False, "retrieval_date": retrieved, "source_release": release, "provenance_reference": "data/processed/xai_reports/dataset_provenance.json",
    }


def _row_decoy(record: dict[str, Any], retrieved: str) -> dict[str, Any]:
    empty = {key: None for key in ("source_assay_id", "source_document_id", "activity_type", "original_activity_value", "original_activity_units", "standard_activity_value_nm", "pchembl_value", "activity_relation", "activity_min_pchembl", "activity_max_pchembl", "activity_median_pchembl")}
    return {
        "compound_id": f"DUDE_{record['source_molecule_id']}", "standardized_compound_id": record["inchikey"], "source_database": "DUD-E", "source_record_id": record["source_record_id"],
        "source_molecule_id": record["source_molecule_id"], "source_target_id": "HIVPR", **empty, "original_smiles": record["original_smiles"],
        **{key: record.get(key) for key in ("canonical_isomeric_smiles", "inchi", "inchikey", "inchikey_connectivity", "molecular_formula", "molecular_weight", "exact_molecular_weight", "logp", "tpsa", "hbd", "hba", "rotatable_bonds", "ring_count", "aromatic_ring_count", "formal_charge", "heavy_atom_count", "fraction_csp3", "qed", "pains_alert_count", "brenk_alert_count", "murcko_scaffold_smiles", "scaffold_id")},
        "binary_label": "inactive_or_decoy", "label_source": "dude_decoy", "activity_record_count": 0, "activity_conflict_flag": False,
        "retrieval_date": retrieved, "source_release": "DUD-E current target files; release not reported", "provenance_reference": "data/processed/xai_reports/dataset_provenance.json",
    }


def main() -> int:
    RDLogger.DisableLog("rdApp.warning")
    args = arguments(); settings = yaml.safe_load((ROOT / "config/settings.yaml").read_text()); config = settings["xai_dataset"]; dirs = settings["directories"]
    paths = {name: ROOT / dirs[name] for name in ("xai_chembl_raw", "xai_dude_raw", "xai_dataset", "xai_reports", "xai_exclusions")}
    for path in paths.values(): path.mkdir(parents=True, exist_ok=True)
    retrieval_metadata_path = ROOT / "data/raw/xai/retrieval_metadata.json"
    if retrieval_metadata_path.is_file(): retrieval_metadata = json.loads(retrieval_metadata_path.read_text(encoding="utf-8"))
    else:
        retrieval_metadata = {"retrieved_at": retrieval_timestamp()}; retrieval_metadata_path.parent.mkdir(parents=True, exist_ok=True); retrieval_metadata_path.write_text(json.dumps(retrieval_metadata, indent=2), encoding="utf-8")
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    chembl = ChEMBLClient(paths["xai_chembl_raw"], config["request_timeout_seconds"], config["request_retries"], config["request_page_size"])
    raw_activities = chembl.fetch_activities(config["target_chembl_id"], args.refresh, args.offline); status = chembl.fetch_status(args.offline)
    release = str(status.get("chembl_db_version") or status.get("version") or "not reported by API")
    curated = curate_activities(raw_activities, config)
    active_outliers = [row for row in curated["active"] if row["molecular_weight"] > config["maximum_molecular_weight"] or row["heavy_atom_count"] < config["minimum_heavy_atoms"] or not config["minimum_formal_charge"] <= row["formal_charge"] <= config["maximum_formal_charge"]]
    active_candidates = [row for row in curated["active"] if row not in active_outliers]
    actives = sorted(active_candidates, key=lambda row: (-row["activity_median_pchembl"], row["inchikey"]))[:args.active_count]
    dude = DUDEClient(paths["xai_dude_raw"], config["request_timeout_seconds"], config["request_retries"])
    dude_files = dude.download(config["dude_target_id"], args.refresh, args.offline); raw_decoys = dude.parse_ism(dude_files["decoys_final.ism"])
    valid_decoys, invalid_decoys, outliers = [], [], []
    for record in raw_decoys:
        try: standardized = {**record, **standardize_structure(record["original_smiles"])}
        except StructureStandardizationError as exc: invalid_decoys.append({**record, "exclusion_reason": str(exc)}); continue
        if standardized["molecular_weight"] > config["maximum_molecular_weight"] or standardized["heavy_atom_count"] < config["minimum_heavy_atoms"] or not config["minimum_formal_charge"] <= standardized["formal_charge"] <= config["maximum_formal_charge"]:
            outliers.append({**standardized, "exclusion_reason": "configured_property_outlier"}); continue
        valid_decoys.append(standardized)
    valid_decoys, duplicate_decoys = resolve_duplicates(valid_decoys, {row["inchikey"] for row in actives})
    selected_decoys, unselected = select_matched_decoys(actives, valid_decoys, min(args.decoy_count, len(actives), args.target_count - len(actives)), args.seed)
    retrieved_at = retrieval_metadata["retrieved_at"]; retrieved = retrieved_at[:10]
    active_rows = [_row_active(row, release, retrieved) for row in actives]; decoy_rows = [_row_decoy(row, retrieved) for row in selected_decoys]
    records = assign_scaffold_splits(active_rows + decoy_rows, args.seed, (config["train_fraction"], config["validation_fraction"], config["test_fraction"]))
    scaffold_frequency = Counter(row["scaffold_id"] for row in records)
    records = [{**row, "scaffold_frequency": scaffold_frequency[row["scaffold_id"]]} for row in records]
    validation = validate_records(records)
    prefix = args.output_prefix; dataset_dir = paths["xai_dataset"]
    csv_path = dataset_dir / f"{prefix}_compounds.csv"; parquet_path = dataset_dir / f"{prefix}_compounds.parquet"; sdf_path = dataset_dir / f"{prefix}_compounds.sdf"; split_path = dataset_dir / "hivpr_xai_scaffold_split.csv"
    frame = pd.DataFrame(records); frame.to_csv(csv_path, index=False); frame.to_parquet(parquet_path, index=False); write_sdf(sdf_path, records)
    write_table(split_path, [{"standardized_compound_id": row["standardized_compound_id"], "scaffold_id": row["scaffold_id"], "dataset_split": row["dataset_split"], "binary_label": row["binary_label"]} for row in records])
    reason_rows = lambda reason: [row for row in curated["invalid"] if row.get("exclusion_reason") == reason]
    exclusions = {"invalid_structures.csv": curated["invalid"] + invalid_decoys, "outside_typical_range.csv": reason_rows("outside_typical_range"), "invalid_measurements.csv": reason_rows("invalid_measurement"), "conflicting_assay_annotations.csv": reason_rows("conflicting_assay_annotation"), "ambiguous_activities.csv": curated["ambiguous"] + reason_rows("ambiguous_measurement_annotation"), "conflicting_activities.csv": curated["conflicting"], "duplicate_compounds.csv": curated["duplicates"] + duplicate_decoys, "property_outliers.csv": active_outliers + outliers, "excluded_decoys.csv": unselected}
    for filename, rows in exclusions.items(): write_table(paths["xai_exclusions"] / filename, rows)
    reports = generate_dataset_reports(records, paths["xai_reports"], active_rows, decoy_rows)
    provenance = write_provenance(paths["xai_reports"] / "dataset_provenance.json", [{"database": "ChEMBL", "target_id": "CHEMBL243", "endpoint": chembl.endpoint, "release": release}, {"database": "DUD-E", "target_id": "HIVPR", "base_url": dude.base_url, "files": [{"filename": name, "sha256": dude.checksum(path), "size_bytes": path.stat().st_size} for name, path in dude_files.items()]}], release, retrieved_at); reports.append(provenance)
    validate_output_files(csv_path, parquet_path, sdf_path)
    stats = {"target_chembl_id": "CHEMBL243", "dude_target_id": "HIVPR", "chembl_release": release, "raw_chembl_activity_count": len(raw_activities), "valid_chembl_activity_record_count": len(curated["valid_records"]), "unique_active_compound_count": len(active_rows), "raw_dude_decoy_count": len(raw_decoys), "valid_dude_decoy_count": len(valid_decoys), "selected_decoy_count": len(decoy_rows), "duplicate_count": len(curated["duplicates"]) + len(duplicate_decoys), "conflicting_activity_count": len(curated["conflicting"]), "exclusion_counts": {name: len(rows) for name, rows in exclusions.items()}, "final_dataset_size": len(records), "class_counts": dict(Counter(row["binary_label"] for row in records)), "unique_scaffold_count": validation["scaffold_count"], "split_counts": dict(Counter(row["dataset_split"] for row in records)), "scaffold_leakage_count": 0}
    generated = [csv_path, parquet_path, sdf_path, split_path] + reports
    manifest, schema, checksums = write_manifest_and_schema(dataset_dir, records, stats, list(dude_files.values()), generated, args.seed)
    print(json.dumps({**stats, "output_paths": [str(path.relative_to(ROOT)) for path in generated + [manifest, schema, checksums]]}, indent=2)); return 0


if __name__ == "__main__": raise SystemExit(main())
