"""Prepare and audit an independent dataset for locked XAI model validation.

This script:
- validates the external dataset schema;
- standardizes molecular structures with RDKit;
- removes invalid and duplicate molecules;
- detects exact overlap with the DockForgeAI training benchmark;
- detects Bemis–Murcko scaffold overlap;
- preserves experimental provenance;
- writes a clean external-validation dataset.

It does not train, tune, calibrate, or modify any model.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold


PROJECT_ROOT = Path(__file__).resolve().parents[1]

DEFAULT_TRAINING_DATASET = (
    PROJECT_ROOT
    / "data/processed/xai_dataset/hivpr_xai_2000_compounds.csv"
)

DEFAULT_OUTPUT_DIRECTORY = (
    PROJECT_ROOT
    / "data/external/xai_validation_prepared"
)


class ExternalDatasetPreparationError(RuntimeError):
    """Raised when an external validation dataset cannot be prepared."""


def canonicalize_smiles(smiles: Any) -> dict[str, Any] | None:
    """Standardize one SMILES string and calculate identity fields."""

    if pd.isna(smiles):
        return None

    value = str(smiles).strip()

    if not value:
        return None

    molecule = Chem.MolFromSmiles(value)

    if molecule is None:
        return None

    fragments = Chem.GetMolFrags(
        molecule,
        asMols=True,
        sanitizeFrags=True,
    )

    organic_fragments = [
        fragment
        for fragment in fragments
        if any(atom.GetAtomicNum() == 6 for atom in fragment.GetAtoms())
    ]

    if organic_fragments:
        molecule = max(
            organic_fragments,
            key=lambda mol: mol.GetNumHeavyAtoms(),
        )
    elif fragments:
        molecule = max(
            fragments,
            key=lambda mol: mol.GetNumHeavyAtoms(),
        )

    Chem.SanitizeMol(molecule)
    Chem.AssignStereochemistry(
        molecule,
        cleanIt=True,
        force=True,
    )

    canonical_smiles = Chem.MolToSmiles(
        molecule,
        canonical=True,
        isomericSmiles=True,
    )

    try:
        inchi = Chem.MolToInchi(molecule)
        inchikey = Chem.InchiToInchiKey(inchi)
    except Exception:
        inchi = ""
        inchikey = ""

    scaffold_molecule = MurckoScaffold.GetScaffoldForMol(molecule)

    scaffold_smiles = (
        Chem.MolToSmiles(
            scaffold_molecule,
            canonical=True,
            isomericSmiles=True,
        )
        if scaffold_molecule.GetNumAtoms() > 0
        else ""
    )

    return {
        "canonical_isomeric_smiles": canonical_smiles,
        "inchi": inchi,
        "inchikey": inchikey,
        "inchikey_connectivity": (
            inchikey.split("-")[0]
            if inchikey
            else ""
        ),
        "murcko_scaffold_smiles": scaffold_smiles,
        "heavy_atom_count": molecule.GetNumHeavyAtoms(),
    }


def normalize_label(value: Any) -> str | None:
    """Map common binary labels to the existing model labels."""

    if pd.isna(value):
        return None

    normalized = str(value).strip().lower()

    positive_values = {
        "1",
        "active",
        "positive",
        "yes",
        "true",
    }

    negative_values = {
        "0",
        "inactive",
        "negative",
        "no",
        "false",
    }

    if normalized in positive_values:
        return "active"

    if normalized in negative_values:
        return "inactive_or_decoy"

    return None


def load_training_identity(
    training_path: Path,
) -> tuple[set[str], set[str], set[str]]:
    """Load training structure and scaffold identities."""

    if not training_path.exists():
        raise FileNotFoundError(
            f"Training dataset not found: {training_path}"
        )

    training = pd.read_csv(training_path)

    required = {
        "canonical_isomeric_smiles",
        "inchikey",
        "murcko_scaffold_smiles",
    }

    missing = required - set(training.columns)

    if missing:
        raise ExternalDatasetPreparationError(
            "Training dataset is missing required columns: "
            f"{sorted(missing)}"
        )

    smiles = set(
        training["canonical_isomeric_smiles"]
        .dropna()
        .astype(str)
    )

    inchikeys = set(
        training["inchikey"]
        .dropna()
        .astype(str)
    ) - {""}

    scaffolds = set(
        training["murcko_scaffold_smiles"]
        .dropna()
        .astype(str)
    ) - {""}

    return smiles, inchikeys, scaffolds


def prepare_external_dataset(
    input_path: Path,
    output_directory: Path,
    training_path: Path,
    id_column: str,
    smiles_column: str,
    label_column: str,
    source_column: str | None,
) -> dict[str, Any]:
    """Prepare an independent external-validation dataset."""

    if not input_path.exists():
        raise FileNotFoundError(
            f"External dataset not found: {input_path}"
        )

    if input_path.suffix.lower() == ".parquet":
        raw = pd.read_parquet(input_path)
    else:
        raw = pd.read_csv(input_path)

    required = {
        id_column,
        smiles_column,
        label_column,
    }

    missing = required - set(raw.columns)

    if missing:
        raise ExternalDatasetPreparationError(
            "External dataset is missing required columns: "
            f"{sorted(missing)}"
        )

    if not output_directory.is_absolute():
        output_directory = (
            PROJECT_ROOT / output_directory
        ).resolve()
    else:
        output_directory = output_directory.resolve()

    output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    training_smiles, training_inchikeys, training_scaffolds = (
        load_training_identity(training_path)
    )

    prepared_rows: list[dict[str, Any]] = []
    excluded_rows: list[dict[str, Any]] = []

    for row_number, row in raw.iterrows():
        compound_id = str(row[id_column]).strip()

        structure = canonicalize_smiles(
            row[smiles_column]
        )

        normalized_label = normalize_label(
            row[label_column]
        )

        base_record = {
            "source_row": int(row_number),
            "compound_id": compound_id,
            "original_smiles": row[smiles_column],
            "original_label": row[label_column],
        }

        if not compound_id:
            excluded_rows.append(
                {
                    **base_record,
                    "exclusion_reason": "missing_compound_id",
                }
            )
            continue

        if structure is None:
            excluded_rows.append(
                {
                    **base_record,
                    "exclusion_reason": "invalid_smiles",
                }
            )
            continue

        if normalized_label is None:
            excluded_rows.append(
                {
                    **base_record,
                    "exclusion_reason": "invalid_binary_label",
                }
            )
            continue

        source_value = (
            str(row[source_column]).strip()
            if source_column and source_column in raw.columns
            else "independent_external_dataset"
        )

        exact_smiles_overlap = (
            structure["canonical_isomeric_smiles"]
            in training_smiles
        )

        exact_inchikey_overlap = (
            bool(structure["inchikey"])
            and structure["inchikey"]
            in training_inchikeys
        )

        scaffold_overlap = (
            bool(structure["murcko_scaffold_smiles"])
            and structure["murcko_scaffold_smiles"]
            in training_scaffolds
        )

        prepared_rows.append(
            {
                "compound_id": compound_id,
                "canonical_smiles": structure[
                    "canonical_isomeric_smiles"
                ],
                "true_label": normalized_label,
                "source_database": source_value,
                "inchi": structure["inchi"],
                "inchikey": structure["inchikey"],
                "inchikey_connectivity": structure[
                    "inchikey_connectivity"
                ],
                "murcko_scaffold_smiles": structure[
                    "murcko_scaffold_smiles"
                ],
                "heavy_atom_count": structure[
                    "heavy_atom_count"
                ],
                "exact_smiles_training_overlap": (
                    exact_smiles_overlap
                ),
                "exact_inchikey_training_overlap": (
                    exact_inchikey_overlap
                ),
                "scaffold_training_overlap": scaffold_overlap,
            }
        )

    prepared = pd.DataFrame(prepared_rows)
    excluded = pd.DataFrame(excluded_rows)

    if prepared.empty:
        raise ExternalDatasetPreparationError(
            "No valid external compounds remained after preparation."
        )

    duplicate_mask = prepared.duplicated(
        subset=["inchikey"],
        keep="first",
    ) & prepared["inchikey"].ne("")

    duplicate_rows = prepared[
        duplicate_mask
    ].copy()

    duplicate_rows["exclusion_reason"] = (
        "duplicate_external_structure"
    )

    prepared = prepared[
        ~duplicate_mask
    ].copy()

    exact_overlap_mask = (
        prepared["exact_smiles_training_overlap"]
        | prepared["exact_inchikey_training_overlap"]
    )

    overlap_rows = prepared[
        exact_overlap_mask
    ].copy()

    overlap_rows["exclusion_reason"] = (
        "exact_training_structure_overlap"
    )

    independent = prepared[
        ~exact_overlap_mask
    ].copy()

    independent = independent.sort_values(
        ["compound_id", "inchikey"],
        kind="stable",
    ).reset_index(drop=True)

    prepared_path = (
        output_directory
        / "external_validation_prepared.csv"
    )

    overlap_path = (
        output_directory
        / "external_training_overlap.csv"
    )

    scaffold_path = (
        output_directory
        / "external_scaffold_overlap.csv"
    )

    excluded_path = (
        output_directory
        / "external_exclusions.csv"
    )

    manifest_path = (
        output_directory
        / "external_preparation_manifest.json"
    )

    independent.to_csv(
        prepared_path,
        index=False,
    )

    overlap_rows.to_csv(
        overlap_path,
        index=False,
    )

    independent[
        independent["scaffold_training_overlap"]
    ].to_csv(
        scaffold_path,
        index=False,
    )

    all_exclusions = pd.concat(
        [
            excluded,
            duplicate_rows,
            overlap_rows,
        ],
        ignore_index=True,
        sort=False,
    )

    all_exclusions.to_csv(
        excluded_path,
        index=False,
    )

    class_counts = (
        independent["true_label"]
        .value_counts()
        .to_dict()
    )

    manifest = {
        "generated_at": datetime.now(
            timezone.utc
        ).isoformat(),
        "input_path": input_path.name,
        "training_dataset": training_path.name,
        "raw_row_count": int(len(raw)),
        "valid_before_overlap_filtering": int(
            len(prepared)
        ),
        "final_independent_row_count": int(
            len(independent)
        ),
        "exact_training_overlap_count": int(
            len(overlap_rows)
        ),
        "external_duplicate_count": int(
            len(duplicate_rows)
        ),
        "invalid_or_excluded_count": int(
            len(excluded)
        ),
        "scaffold_overlap_count": int(
            independent[
                "scaffold_training_overlap"
            ].astype(bool).sum()
        ),
        "novel_scaffold_count": int(
            (~independent[
                "scaffold_training_overlap"
            ]).sum()
        ),
        "class_counts": class_counts,
        "model_retraining_performed": False,
        "threshold_tuning_performed": False,
        "output_files": {
            "prepared_dataset": str(
                prepared_path.relative_to(
                    PROJECT_ROOT
                )
            ),
            "training_overlap": str(
                overlap_path.relative_to(
                    PROJECT_ROOT
                )
            ),
            "scaffold_overlap": str(
                scaffold_path.relative_to(
                    PROJECT_ROOT
                )
            ),
            "exclusions": str(
                excluded_path.relative_to(
                    PROJECT_ROOT
                )
            ),
        },
        "limitations": [
            (
                "This preparation step does not establish "
                "biological validity."
            ),
            (
                "Scaffold-overlapping compounds are retained and "
                "flagged, but exact training structures are removed."
            ),
            (
                "The locked model must not be retrained or tuned "
                "using this external dataset."
            ),
        ],
    }

    manifest_path.write_text(
        json.dumps(
            manifest,
            indent=2,
        ),
        encoding="utf-8",
    )

    return manifest


def build_parser() -> argparse.ArgumentParser:
    """Create command-line parser."""

    parser = argparse.ArgumentParser(
        description=(
            "Prepare an independent molecular dataset for "
            "locked-model external validation."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--training-dataset",
        type=Path,
        default=DEFAULT_TRAINING_DATASET,
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIRECTORY,
    )
    parser.add_argument(
        "--id-column",
        default="compound_id",
    )
    parser.add_argument(
        "--smiles-column",
        default="canonical_smiles",
    )
    parser.add_argument(
        "--label-column",
        default="true_label",
    )
    parser.add_argument(
        "--source-column",
        default="source_database",
    )

    return parser


def main() -> int:
    """Run external-dataset preparation."""

    parser = build_parser()
    arguments = parser.parse_args()

    try:
        manifest = prepare_external_dataset(
            input_path=arguments.input,
            output_directory=arguments.output_dir,
            training_path=arguments.training_dataset,
            id_column=arguments.id_column,
            smiles_column=arguments.smiles_column,
            label_column=arguments.label_column,
            source_column=arguments.source_column,
        )
    except (
        FileNotFoundError,
        ExternalDatasetPreparationError,
        ValueError,
    ) as error:
        print(f"External preparation failed: {error}")
        return 1

    print("External dataset preparation completed.")
    print(
        "Final independent compounds:",
        manifest["final_independent_row_count"],
    )
    print(
        "Exact training overlaps removed:",
        manifest["exact_training_overlap_count"],
    )
    print(
        "Scaffold-overlapping compounds:",
        manifest["scaffold_overlap_count"],
    )
    print(
        "Novel-scaffold compounds:",
        manifest["novel_scaffold_count"],
    )
    print(
        "Class counts:",
        manifest["class_counts"],
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
