"""Summary tables, schema, manifest, SDF, and checksums for XAI-1."""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem

from src.datasets.decoy_matcher import matching_report
from src.datasets.provenance import sha256_file


def write_table(path: Path, rows: list[dict[str, Any]], columns: list[str] | None = None) -> Path:
    """Write a CSV with stable columns, including for empty exclusions."""
    path.parent.mkdir(parents=True, exist_ok=True); pd.DataFrame(rows, columns=columns).to_csv(path, index=False); return path


def write_sdf(path: Path, records: list[dict[str, Any]]) -> Path:
    """Write final standardized molecules with dataset metadata."""
    path.parent.mkdir(parents=True, exist_ok=True); writer = Chem.SDWriter(str(path))
    try:
        for row in records:
            molecule = Chem.MolFromSmiles(row["canonical_isomeric_smiles"])
            if molecule is None: continue
            for key in ("compound_id", "standardized_compound_id", "binary_label", "label_source", "dataset_split", "scaffold_id"): molecule.SetProp(key, str(row[key]))
            writer.write(molecule)
    finally: writer.close()
    return path


def generate_dataset_reports(records: list[dict[str, Any]], report_dir: Path, actives: list[dict[str, Any]], decoys: list[dict[str, Any]]) -> list[Path]:
    """Generate deterministic dataset distribution and matching reports."""
    report_dir.mkdir(parents=True, exist_ok=True); outputs = []
    def counts(column: str, filename: str) -> None:
        outputs.append(write_table(report_dir / filename, [{column: key, "count": value} for key, value in sorted(Counter(str(row.get(column, "")) for row in records).items())]))
    counts("binary_label", "class_distribution.csv"); counts("source_database", "source_distribution.csv"); counts("scaffold_id", "scaffold_distribution.csv"); counts("dataset_split", "split_distribution.csv")
    outputs.append(write_table(report_dir / "activity_distribution.csv", [{"compound_id": row["compound_id"], "pchembl_value": row.get("pchembl_value"), "standard_activity_value_nm": row.get("standard_activity_value_nm")} for row in records if row["binary_label"] == "active"]))
    descriptor_columns = ("molecular_weight", "logp", "tpsa", "hbd", "hba", "rotatable_bonds", "formal_charge")
    descriptor_rows = []
    frame = pd.DataFrame(records)
    for label, group in frame.groupby("binary_label"):
        for column in descriptor_columns: descriptor_rows.append({"binary_label": label, "descriptor": column, "mean": group[column].mean(), "median": group[column].median(), "minimum": group[column].min(), "maximum": group[column].max()})
    outputs.append(write_table(report_dir / "descriptor_distribution.csv", descriptor_rows))
    outputs.append(write_table(report_dir / "active_decoy_matching.csv", matching_report(actives, decoys)))
    summary = {"dataset_size": len(records), "class_counts": Counter(row["binary_label"] for row in records), "source_counts": Counter(row["source_database"] for row in records), "unique_scaffolds": len({row["scaffold_id"] for row in records}), "split_counts": Counter(row["dataset_split"] for row in records), "scientific_limitations": ["DUD-E decoys are computational property-matched decoys, not experimentally proven inactive compounds.", "This dataset phase does not include docking, feature engineering, machine learning, or model performance claims."]}
    summary_json = report_dir / "dataset_summary.json"; summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8"); outputs.append(summary_json)
    summary_md = report_dir / "dataset_summary.md"; summary_md.write_text(f"# HIVPR XAI-1 dataset summary\n\n- Final compounds: {len(records)}\n- Active compounds: {summary['class_counts'].get('active', 0)}\n- DUD-E decoys: {summary['class_counts'].get('inactive_or_decoy', 0)}\n- Unique scaffolds: {summary['unique_scaffolds']}\n\nDUD-E decoys are computational decoys, not experimentally proven inactive compounds. No model was trained.\n", encoding="utf-8"); outputs.append(summary_md)
    return outputs


def write_manifest_and_schema(dataset_dir: Path, records: list[dict[str, Any]], stats: dict[str, Any], source_files: list[Path], generated_files: list[Path], seed: int) -> tuple[Path, Path, Path]:
    """Write manifest, schema, and checksums for final artifacts."""
    dataset_dir.mkdir(parents=True, exist_ok=True)
    schema_path = dataset_dir / "hivpr_xai_dataset_schema.json"
    schema = {"columns": [{"name": column, "dtype": str(pd.DataFrame(records)[column].dtype), "nullable": bool(pd.DataFrame(records)[column].isna().any())} for column in pd.DataFrame(records).columns], "binary_label_values": ["active", "inactive_or_decoy"], "dataset_split_values": ["train", "validation", "test"]}
    schema_path.write_text(json.dumps(schema, indent=2), encoding="utf-8")
    checksums_path = dataset_dir / "hivpr_xai_dataset_checksums.json"
    checksum_files = [path for path in generated_files + [schema_path] if path.is_file()]
    checksums_path.write_text(json.dumps({path.name: {"sha256": sha256_file(path), "size_bytes": path.stat().st_size} for path in checksum_files}, indent=2), encoding="utf-8")
    manifest_path = dataset_dir / "hivpr_xai_dataset_manifest.json"
    manifest = {**stats, "random_seed": seed, "deterministic_with_cached_inputs": True, "source_files": [{"filename": path.name, "sha256": sha256_file(path), "size_bytes": path.stat().st_size} for path in source_files if path.is_file()], "generated_files": [path.name for path in generated_files], "schema_file": schema_path.name, "checksums_file": checksums_path.name}
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest_path, schema_path, checksums_path
