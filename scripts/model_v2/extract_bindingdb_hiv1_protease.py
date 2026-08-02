"""Stream-filter the complete BindingDB TSV for HIV-1 protease records.

This script does not train a model or assign activity labels. It creates a
small, target-specific CSV with provenance and measurement columns retained.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Iterable

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

TARGET_PATTERN = re.compile(
    r"("
    r"\bHIV[\s_-]*1\b.*\bprotease\b"
    r"|human immunodeficiency virus(?:\s+type\s+1)?\s+protease"
    r"|\bHIV\s+protease\b"
    r")",
    flags=re.IGNORECASE,
)

MUTATION_PATTERN = re.compile(
    r"("
    r"\bmutant\b|\bvariant\b|\bmutation\b"
    r"|\b[A-Z]\d{1,4}[A-Z]\b"
    r")",
    flags=re.IGNORECASE,
)


def select_columns(columns: Iterable[str]) -> list[str]:
    """Select target, ligand, measurement and provenance columns."""

    selected: list[str] = []

    important_terms = (
        "smiles",
        "inchi",
        "monomerid",
        "reactant_set_id",
        "ligand",
        "target",
        "uniprot",
        "pdb",
        "sequence",
        "organism",
        "ki",
        "kd",
        "ic50",
        "ec50",
        "kon",
        "koff",
        "doi",
        "pubmed",
        "patent",
        "article",
        "journal",
        "authors",
        "institution",
        "assay",
        "measurement",
        "entry name",
        "primary id",
        "secondary id",
        "alternative id",
    )

    for column in columns:
        lowered = column.lower()

        if any(term in lowered for term in important_terms):
            selected.append(column)

    return list(dict.fromkeys(selected))


def target_search_columns(columns: Iterable[str]) -> list[str]:
    """Find columns likely to describe target identity."""

    output: list[str] = []

    for column in columns:
        lowered = column.lower()

        if (
            "target" in lowered
            or "uniprot" in lowered
            or "organism" in lowered
        ):
            output.append(column)

    return output


def normalize_text_frame(frame: pd.DataFrame, columns: list[str]) -> pd.Series:
    """Concatenate target metadata safely for text matching."""

    available = [column for column in columns if column in frame.columns]

    if not available:
        return pd.Series("", index=frame.index, dtype="string")

    return (
        frame[available]
        .fillna("")
        .astype(str)
        .agg(" | ".join, axis=1)
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Extract HIV-1 protease rows from a complete BindingDB TSV."
    )
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "data/model_v2/raw/"
            "bindingdb_hiv1_protease_filtered.csv"
        ),
    )
    parser.add_argument(
        "--audit-output",
        type=Path,
        default=Path(
            "data/model_v2/audit/"
            "bindingdb_hiv1_protease_extraction.json"
        ),
    )
    parser.add_argument("--chunksize", type=int, default=2500)
    args = parser.parse_args()

    input_path = args.input
    output_path = args.output
    audit_path = args.audit_output

    if not input_path.is_absolute():
        input_path = (PROJECT_ROOT / input_path).resolve()

    if not output_path.is_absolute():
        output_path = (PROJECT_ROOT / output_path).resolve()

    if not audit_path.is_absolute():
        audit_path = (PROJECT_ROOT / audit_path).resolve()

    if not input_path.exists():
        raise FileNotFoundError(f"BindingDB TSV not found: {input_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.parent.mkdir(parents=True, exist_ok=True)

    # Read only the header first.
    header = pd.read_csv(
        input_path,
        sep="\t",
        nrows=0,
        dtype=str,
    )

    all_columns = list(header.columns)
    selected_columns = select_columns(all_columns)
    search_columns = target_search_columns(selected_columns)

    if not selected_columns:
        raise RuntimeError("No usable BindingDB columns were detected.")

    if not search_columns:
        raise RuntimeError("No target-identity columns were detected.")

    print("Input:", input_path)
    print("Total BindingDB columns:", len(all_columns))
    print("Selected columns:", len(selected_columns))
    print("Target-search columns:", len(search_columns))
    print("Chunk size:", args.chunksize)

    if output_path.exists():
        output_path.unlink()

    total_rows = 0
    matched_rows = 0
    possible_mutant_rows = 0
    chunk_number = 0
    wrote_header = False

    reader = pd.read_csv(
        input_path,
        sep="\t",
        dtype=str,
        usecols=selected_columns,
        chunksize=args.chunksize,
        low_memory=False,
        on_bad_lines="warn",
        encoding="utf-8",
        encoding_errors="replace",
    )

    for chunk in reader:
        chunk_number += 1
        total_rows += len(chunk)

        target_text = normalize_text_frame(chunk, search_columns)
        match_mask = target_text.str.contains(
            TARGET_PATTERN,
            na=False,
        )

        matched = chunk.loc[match_mask].copy()

        if not matched.empty:
            matched["_bindingdb_target_text"] = target_text.loc[
                match_mask
            ].values

            matched["_possible_mutant_or_variant"] = (
                matched["_bindingdb_target_text"]
                .str.contains(MUTATION_PATTERN, na=False)
            )

            possible_mutant_rows += int(
                matched["_possible_mutant_or_variant"].sum()
            )
            matched_rows += len(matched)

            matched.to_csv(
                output_path,
                mode="a",
                header=not wrote_header,
                index=False,
            )
            wrote_header = True

        if chunk_number % 25 == 0:
            print(
                f"Chunks: {chunk_number:,} | "
                f"Rows scanned: {total_rows:,} | "
                f"Matches: {matched_rows:,}",
                flush=True,
            )

    if not wrote_header:
        pd.DataFrame(columns=selected_columns).to_csv(
            output_path,
            index=False,
        )

    audit = {
        "input_file": input_path.name,
        "input_size_bytes": input_path.stat().st_size,
        "total_columns": len(all_columns),
        "selected_column_count": len(selected_columns),
        "target_search_column_count": len(search_columns),
        "chunks_processed": chunk_number,
        "rows_scanned": total_rows,
        "matched_rows": matched_rows,
        "possible_mutant_or_variant_rows": possible_mutant_rows,
        "output_file": str(output_path.relative_to(PROJECT_ROOT)),
        "selected_columns": selected_columns,
        "target_search_columns": search_columns,
        "notes": [
            "Records are target-text matches only.",
            "Possible mutants and variants are flagged, not silently removed.",
            "No activity labels were assigned.",
            "No model training was performed.",
        ],
    }

    audit_path.write_text(
        json.dumps(audit, indent=2),
        encoding="utf-8",
    )

    print("\nExtraction complete")
    print("Rows scanned:", total_rows)
    print("Matched rows:", matched_rows)
    print("Possible mutant/variant rows:", possible_mutant_rows)
    print("Output:", output_path)
    print("Audit:", audit_path)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
