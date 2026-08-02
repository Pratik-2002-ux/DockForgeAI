"""Audit the filtered BindingDB HIV-1 protease dataset before curation."""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem


PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT = (
    PROJECT_ROOT
    / "data/model_v2/raw/bindingdb_hiv1_protease_filtered.csv"
)

OUTPUT_DIR = PROJECT_ROOT / "data/model_v2/audit"


def parse_measurement(value: object) -> tuple[str | None, float | None]:
    """Parse BindingDB values such as 0.31, <0.01, >10000."""

    if pd.isna(value):
        return None, None

    text = str(value).strip()

    if not text:
        return None, None

    match = re.fullmatch(
        r"\s*(<=|>=|<|>|=)?\s*"
        r"([0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)\s*",
        text,
    )

    if not match:
        return None, None

    relation = match.group(1) or "="
    numeric = float(match.group(2))

    return relation, numeric


def canonicalize(smiles: object) -> str | None:
    """Return canonical isomeric SMILES or None."""

    if pd.isna(smiles):
        return None

    mol = Chem.MolFromSmiles(str(smiles))

    if mol is None:
        return None

    fragments = Chem.GetMolFrags(
        mol,
        asMols=True,
        sanitizeFrags=True,
    )

    if fragments:
        organic = [
            fragment
            for fragment in fragments
            if any(
                atom.GetAtomicNum() == 6
                for atom in fragment.GetAtoms()
            )
        ]

        mol = max(
            organic or list(fragments),
            key=lambda item: item.GetNumHeavyAtoms(),
        )

    return Chem.MolToSmiles(
        mol,
        canonical=True,
        isomericSmiles=True,
    )


def main() -> int:
    if not INPUT.exists():
        raise FileNotFoundError(INPUT)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(
        INPUT,
        dtype=str,
        low_memory=False,
    )

    measurement_columns = [
        column
        for column in [
            "Ki (nM)",
            "Kd (nM)",
            "IC50 (nM)",
            "EC50 (nM)",
        ]
        if column in df.columns
    ]

    records: list[dict[str, object]] = []

    for row_index, row in df.iterrows():
        canonical = canonicalize(row.get("Ligand SMILES"))

        for activity_type in measurement_columns:
            relation, value_nm = parse_measurement(
                row.get(activity_type)
            )

            if value_nm is None:
                continue

            if relation in {"=", "<", "<="} and value_nm <= 1000:
                label = "active"
            elif relation in {"=", ">", ">="} and value_nm >= 10000:
                label = "inactive"
            else:
                label = "ambiguous"

            records.append(
                {
                    "source_row": int(row_index),
                    "bindingdb_reactant_set_id": row.get(
                        "BindingDB Reactant_set_id"
                    ),
                    "bindingdb_monomer_id": row.get(
                        "BindingDB MonomerID"
                    ),
                    "ligand_smiles": row.get("Ligand SMILES"),
                    "canonical_smiles": canonical,
                    "activity_type": activity_type.replace(
                        " (nM)", ""
                    ),
                    "relation": relation,
                    "value_nm": value_nm,
                    "binary_status": label,
                    "article_doi": row.get("Article DOI"),
                    "institution": row.get("Institution"),
                }
            )

    measurements = pd.DataFrame(records)

    measurements.to_csv(
        OUTPUT_DIR
        / "bindingdb_hiv1_measurement_audit.csv",
        index=False,
    )

    valid_structures = measurements[
        measurements["canonical_smiles"].notna()
    ].copy()

    compound_summary = (
        valid_structures.groupby("canonical_smiles")
        .agg(
            measurement_count=("value_nm", "size"),
            min_value_nm=("value_nm", "min"),
            median_value_nm=("value_nm", "median"),
            max_value_nm=("value_nm", "max"),
            active_measurements=(
                "binary_status",
                lambda values: int((values == "active").sum()),
            ),
            inactive_measurements=(
                "binary_status",
                lambda values: int((values == "inactive").sum()),
            ),
            ambiguous_measurements=(
                "binary_status",
                lambda values: int((values == "ambiguous").sum()),
            ),
            activity_types=(
                "activity_type",
                lambda values: ";".join(
                    sorted(set(values.dropna()))
                ),
            ),
            source_dois=(
                "article_doi",
                lambda values: ";".join(
                    sorted(
                        {
                            str(value)
                            for value in values.dropna()
                            if str(value).strip()
                        }
                    )
                ),
            ),
        )
        .reset_index()
    )

    compound_summary["conflict_flag"] = (
        compound_summary["active_measurements"].gt(0)
        & compound_summary["inactive_measurements"].gt(0)
    )

    compound_summary["final_status"] = np.select(
        [
            compound_summary["conflict_flag"],
            (
                compound_summary["active_measurements"].gt(0)
                & compound_summary["inactive_measurements"].eq(0)
            ),
            (
                compound_summary["inactive_measurements"].gt(0)
                & compound_summary["active_measurements"].eq(0)
            ),
        ],
        [
            "conflicting",
            "active",
            "inactive",
        ],
        default="ambiguous_only",
    )

    compound_summary.to_csv(
        OUTPUT_DIR
        / "bindingdb_hiv1_compound_audit.csv",
        index=False,
    )

    summary = {
        "raw_rows": int(len(df)),
        "parsed_measurements": int(len(measurements)),
        "invalid_or_missing_structure_measurements": int(
            measurements["canonical_smiles"].isna().sum()
        ),
        "unique_valid_compounds": int(
            compound_summary["canonical_smiles"].nunique()
        ),
        "measurement_type_counts": {
            str(key): int(value)
            for key, value in (
                measurements["activity_type"]
                .value_counts()
                .to_dict()
                .items()
            )
        },
        "measurement_label_counts": {
            str(key): int(value)
            for key, value in (
                measurements["binary_status"]
                .value_counts()
                .to_dict()
                .items()
            )
        },
        "compound_status_counts": {
            str(key): int(value)
            for key, value in (
                compound_summary["final_status"]
                .value_counts()
                .to_dict()
                .items()
            )
        },
        "unique_doi_count": int(
            measurements["article_doi"]
            .dropna()
            .astype(str)
            .str.strip()
            .replace("", np.nan)
            .dropna()
            .nunique()
        ),
        "training_readiness_rule": {
            "minimum_active_compounds": 200,
            "minimum_inactive_compounds": 200,
            "minimum_total_compounds_for_gnn": 500,
        },
    }

    active_count = summary["compound_status_counts"].get(
        "active", 0
    )
    inactive_count = summary["compound_status_counts"].get(
        "inactive", 0
    )

    summary["preliminary_training_ready"] = bool(
        active_count >= 200
        and inactive_count >= 200
    )

    summary["preliminary_gnn_ready"] = bool(
        summary["unique_valid_compounds"] >= 500
        and summary["preliminary_training_ready"]
    )

    summary_path = (
        OUTPUT_DIR
        / "bindingdb_hiv1_audit_summary.json"
    )

    summary_path.write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(summary, indent=2))
    print("\nMeasurement audit:")
    print(
        OUTPUT_DIR
        / "bindingdb_hiv1_measurement_audit.csv"
    )
    print("Compound audit:")
    print(
        OUTPUT_DIR
        / "bindingdb_hiv1_compound_audit.csv"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
