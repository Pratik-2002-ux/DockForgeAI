"""Build a transparent Phase 4 consensus ranking.

Combines frozen Track 1 predictions, applicability-domain evidence,
Vina docking, pose analysis, interaction evidence, structural alerts,
chemical diversity, and supplier availability.

This script prioritizes candidates for expert review only.
It does not establish experimental inhibition or binding.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold


PROJECT_ROOT = Path(__file__).resolve().parents[1]

DEFAULT_WEIGHTS = {
    "predicted_pki_component": 0.25,
    "docking_component": 0.20,
    "applicability_component": 0.15,
    "pose_quality_component": 0.15,
    "interaction_component": 0.10,
    "uncertainty_component": 0.05,
    "structural_alert_component": 0.05,
    "diversity_component": 0.03,
    "supplier_component": 0.02,
}

SCENARIOS = {
    "balanced_default": DEFAULT_WEIGHTS,
    "ml_emphasis": {
        "predicted_pki_component": 0.40,
        "docking_component": 0.15,
        "applicability_component": 0.15,
        "pose_quality_component": 0.10,
        "interaction_component": 0.08,
        "uncertainty_component": 0.05,
        "structural_alert_component": 0.03,
        "diversity_component": 0.02,
        "supplier_component": 0.02,
    },
    "docking_emphasis": {
        "predicted_pki_component": 0.15,
        "docking_component": 0.40,
        "applicability_component": 0.10,
        "pose_quality_component": 0.15,
        "interaction_component": 0.10,
        "uncertainty_component": 0.03,
        "structural_alert_component": 0.03,
        "diversity_component": 0.02,
        "supplier_component": 0.02,
    },
    "domain_quality_emphasis": {
        "predicted_pki_component": 0.15,
        "docking_component": 0.15,
        "applicability_component": 0.25,
        "pose_quality_component": 0.20,
        "interaction_component": 0.10,
        "uncertainty_component": 0.05,
        "structural_alert_component": 0.05,
        "diversity_component": 0.03,
        "supplier_component": 0.02,
    },
    "equal_weight": {
        key: 1.0 / len(DEFAULT_WEIGHTS)
        for key in DEFAULT_WEIGHTS
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_column(
    frame: pd.DataFrame,
    aliases: Iterable[str],
    required: bool = True,
) -> str | None:
    lowered = {column.lower(): column for column in frame.columns}

    for alias in aliases:
        if alias.lower() in lowered:
            return lowered[alias.lower()]

    if required:
        raise KeyError(
            f"None of these columns were found: {list(aliases)}. "
            f"Available columns: {frame.columns.tolist()}"
        )

    return None


def robust_scale(
    series: pd.Series,
    higher_is_better: bool = True,
) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    valid = values.dropna()

    if valid.empty:
        return pd.Series(np.nan, index=series.index)

    lower = valid.quantile(0.05)
    upper = valid.quantile(0.95)

    if np.isclose(lower, upper):
        scaled = pd.Series(0.5, index=series.index)
        scaled[values.isna()] = np.nan
        return scaled

    clipped = values.clip(lower=lower, upper=upper)
    scaled = (clipped - lower) / (upper - lower)

    if not higher_is_better:
        scaled = 1.0 - scaled

    return scaled.clip(0.0, 1.0)


def as_bool(value: object) -> bool:
    if pd.isna(value):
        return False

    if isinstance(value, bool):
        return value

    return str(value).strip().lower() in {
        "true",
        "1",
        "yes",
        "y",
    }


def availability_score(value: object) -> float:
    text = str(value).strip().lower()

    if "in_stock" in text or "listed_in_stock" in text:
        return 1.0

    if text in {"unknown", "", "nan"}:
        return 0.5

    if "unavailable" in text or "out_of_stock" in text:
        return 0.0

    return 0.5


def applicability_score(value: object) -> float:
    text = str(value).strip().lower()

    if text == "in_domain":
        return 1.0

    if text == "borderline":
        return 0.5

    if text == "out_of_domain":
        return 0.0

    return np.nan


def alert_score(pains: object, brenk: object) -> float:
    pains_flag = as_bool(pains)
    brenk_flag = as_bool(brenk)

    if not pains_flag and not brenk_flag:
        return 1.0

    if brenk_flag and not pains_flag:
        return 0.6

    if pains_flag and not brenk_flag:
        return 0.4

    return 0.2


def pose_quality_score(row: pd.Series) -> float:
    ligand_atoms = float(row.get("ligand_atom_count", 0) or 0)
    contacting_residues = float(
        row.get("contacting_residue_count", 0) or 0
    )

    if ligand_atoms <= 0:
        return 0.0

    contact_density = contacting_residues / max(ligand_atoms, 1.0)

    return float(np.clip(contact_density / 0.35, 0.0, 1.0))


def interaction_score(row: pd.Series) -> float:
    hbond = float(row.get("hydrogen_bond_count", 0) or 0)
    hydrophobic = float(
        row.get("hydrophobic_contact_count", 0) or 0
    )
    salt = float(row.get("salt_bridge_count", 0) or 0)
    pi_stack = float(row.get("pi_stacking_count", 0) or 0)
    pi_cation = float(row.get("pi_cation_count", 0) or 0)
    metal = float(row.get("metal_interaction_count", 0) or 0)

    score = (
        min(hbond, 4.0) / 4.0 * 0.30
        + min(hydrophobic, 12.0) / 12.0 * 0.30
        + min(salt, 2.0) / 2.0 * 0.15
        + min(pi_stack, 2.0) / 2.0 * 0.10
        + min(pi_cation, 2.0) / 2.0 * 0.10
        + min(metal, 1.0) * 0.05
    )

    return float(np.clip(score, 0.0, 1.0))


def fingerprint(smiles: str):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None, ""

    fp = AllChem.GetMorganGenerator(
        radius=2,
        fpSize=2048,
    ).GetFingerprint(mol)

    scaffold_mol = MurckoScaffold.GetScaffoldForMol(mol)
    scaffold = (
        Chem.MolToSmiles(
            scaffold_mol,
            canonical=True,
            isomericSmiles=True,
        )
        if scaffold_mol is not None
        and scaffold_mol.GetNumAtoms() > 0
        else ""
    )

    return fp, scaffold


def calculate_diversity(
    frame: pd.DataFrame,
    smiles_column: str,
) -> pd.DataFrame:
    """Calculate within-pool diversity from one unambiguous SMILES column."""

    if smiles_column not in frame.columns:
        raise KeyError(
            f"SMILES column not found: {smiles_column}"
        )

    selected = frame.loc[:, smiles_column]

    if isinstance(selected, pd.DataFrame):
        raise ValueError(
            f"Duplicate SMILES column detected: {smiles_column}"
        )

    fingerprints = []
    scaffolds = []

    for row_index, smiles in selected.items():
        text_smiles = str(smiles).strip()

        fp, scaffold = fingerprint(text_smiles)

        if fp is None:
            compound_id = frame.loc[row_index, "compound_id"]
            raise ValueError(
                "RDKit could not parse the prospective structure for "
                f"{compound_id}: {text_smiles}"
            )

        fingerprints.append(fp)
        scaffolds.append(scaffold)

    if len(scaffolds) != len(frame):
        raise RuntimeError(
            "Diversity scaffold count does not match candidate count: "
            f"{len(scaffolds)} versus {len(frame)}"
        )

    frame = frame.copy()
    frame["murcko_scaffold_smiles_consensus"] = scaffolds

    similarities = []

    for index, fp in enumerate(fingerprints):
        if fp is None:
            similarities.append(np.nan)
            continue

        other_fps = [
            candidate
            for other_index, candidate in enumerate(fingerprints)
            if other_index != index and candidate is not None
        ]

        if not other_fps:
            similarities.append(0.0)
            continue

        sims = DataStructs.BulkTanimotoSimilarity(fp, other_fps)
        similarities.append(max(sims))

    frame["nearest_candidate_similarity"] = similarities
    frame["diversity_component"] = (
        1.0 - frame["nearest_candidate_similarity"]
    ).clip(0.0, 1.0)

    return frame


def weighted_score(
    frame: pd.DataFrame,
    weights: dict[str, float],
) -> tuple[pd.Series, pd.DataFrame]:
    scores = []
    effective_rows = []

    for _, row in frame.iterrows():
        available = {
            component: weight
            for component, weight in weights.items()
            if pd.notna(row.get(component))
        }

        total_weight = sum(available.values())

        if total_weight <= 0:
            scores.append(np.nan)
            effective_rows.append({})
            continue

        effective = {
            component: weight / total_weight
            for component, weight in available.items()
        }

        score = sum(
            float(row[component]) * weight
            for component, weight in effective.items()
        )

        scores.append(float(score))
        effective_rows.append(effective)

    effective_frame = pd.DataFrame(effective_rows)

    return pd.Series(scores, index=frame.index), effective_frame


def git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=PROJECT_ROOT,
            text=True,
        ).strip()
    except Exception:
        return "unavailable"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--study-id",
        default="HIVPR_2026_001",
    )
    parser.add_argument(
        "--screening",
        type=Path,
        default=Path(
            "data/processed/prospective_studies/"
            "HIVPR_2026_001/ml_screening/"
            "top_30_diverse_review_pool.csv"
        ),
    )
    parser.add_argument(
        "--docking",
        type=Path,
        default=Path(
            "data/processed/phase4/docking/"
            "docking_results.csv"
        ),
    )
    parser.add_argument(
        "--pose-summary",
        type=Path,
        default=Path(
            "data/processed/phase4/analysis/"
            "pose_summary.csv"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "data/processed/phase4/consensus"
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
    )
    args = parser.parse_args()

    paths = [
        args.screening,
        args.docking,
        args.pose_summary,
    ]

    for path in paths:
        if not path.exists():
            raise FileNotFoundError(path)

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    if any(output_dir.iterdir()) and not args.overwrite:
        raise FileExistsError(
            f"{output_dir} is not empty. Use --overwrite."
        )

    screening = pd.read_csv(
        args.screening,
        low_memory=False,
    )
    docking = pd.read_csv(
        args.docking,
        low_memory=False,
    )
    pose = pd.read_csv(
        args.pose_summary,
        low_memory=False,
    )

    screening_id = resolve_column(
        screening,
        ["compound_id"],
    )
    docking_id = resolve_column(
        docking,
        ["compound_id"],
    )
    pose_id = resolve_column(
        pose,
        ["compound_id"],
    )

    screening = screening.rename(
        columns={screening_id: "compound_id"}
    )
    docking = docking.rename(
        columns={docking_id: "compound_id"}
    )
    pose = pose.rename(
        columns={pose_id: "compound_id"}
    )

    score_column = resolve_column(
        docking,
        [
            "docking_score_kcal_mol",
            "vina_score_kcal_mol",
            "best_affinity_kcal_mol",
            "best_score",
            "affinity",
        ],
    )
    docking = docking.rename(
        columns={score_column: "vina_score_kcal_mol"}
    )

    # The docking table contains up to 10 Vina poses per compound.
    # Consensus ranking uses the predefined best pose, pose_rank == 1.
    if "pose_rank" not in docking.columns:
        raise KeyError(
            "docking_results.csv must contain pose_rank so the best pose "
            "can be selected reproducibly."
        )

    docking["pose_rank"] = pd.to_numeric(
        docking["pose_rank"],
        errors="coerce",
    )
    docking["vina_score_kcal_mol"] = pd.to_numeric(
        docking["vina_score_kcal_mol"],
        errors="coerce",
    )

    all_pose_counts = (
        docking.groupby("compound_id")
        .size()
        .rename("docking_pose_count")
        .reset_index()
    )

    docking = docking[
        docking["pose_rank"].eq(1)
        & docking["vina_score_kcal_mol"].notna()
    ].copy()

    docking = docking.merge(
        all_pose_counts,
        on="compound_id",
        how="left",
    )

    if docking["compound_id"].duplicated().any():
        duplicated = (
            docking.loc[
                docking["compound_id"].duplicated(keep=False),
                "compound_id",
            ]
            .astype(str)
            .unique()
            .tolist()
        )
        raise ValueError(
            "Multiple pose_rank=1 records found for compounds: "
            f"{duplicated[:20]}"
        )

    print(
        "Docking table:",
        int(all_pose_counts["docking_pose_count"].sum()),
        "poses across",
        len(all_pose_counts),
        "compounds; selected",
        len(docking),
        "pose-rank-1 records",
    )

    predicted_pki_column = resolve_column(
        screening,
        ["predicted_pki"],
    )
    smiles_column = resolve_column(
        screening,
        [
            "canonical_isomeric_smiles",
            "canonical_smiles",
        ],
    )

    # Copy the selected structure field into a dedicated consensus column.
    # Do not rename it to canonical_smiles because the screening table may
    # already contain both canonical_smiles and canonical_isomeric_smiles,
    # which would create duplicate DataFrame column names.
    screening = screening.copy()

    if predicted_pki_column != "predicted_pki":
        screening["predicted_pki"] = screening[predicted_pki_column]

    selected_smiles = screening[smiles_column]

    if isinstance(selected_smiles, pd.DataFrame):
        raise ValueError(
            f"Selected SMILES column is duplicated: {smiles_column}. "
            "Input columns must be unique."
        )

    screening["consensus_smiles"] = (
        selected_smiles
        .fillna("")
        .astype(str)
        .str.strip()
    )

    invalid_smiles_text = screening["consensus_smiles"].isin(
        ["", "nan", "None", "canonical_smiles"]
    )
    if invalid_smiles_text.any():
        invalid_ids = screening.loc[
            invalid_smiles_text,
            "compound_id",
        ].astype(str).tolist()

        raise ValueError(
            "Missing or malformed screening SMILES for compounds: "
            f"{invalid_ids[:20]}"
        )

    merged = screening.merge(
        docking,
        on="compound_id",
        how="left",
        suffixes=("", "_docking"),
    ).merge(
        pose,
        on="compound_id",
        how="left",
        suffixes=("", "_pose"),
    )

    if merged["compound_id"].duplicated().any():
        duplicates = merged.loc[
            merged["compound_id"].duplicated(False),
            "compound_id",
        ].tolist()
        raise ValueError(
            f"Duplicate candidates after merging: {duplicates}"
        )

    merged["predicted_pki_component"] = robust_scale(
        merged["predicted_pki"],
        higher_is_better=True,
    )
    merged["docking_component"] = robust_scale(
        merged["vina_score_kcal_mol"],
        higher_is_better=False,
    )

    domain_column = resolve_column(
        merged,
        ["applicability_domain_status"],
        required=False,
    )
    if domain_column:
        merged["applicability_component"] = merged[
            domain_column
        ].apply(applicability_score)
    else:
        merged["applicability_component"] = np.nan

    uncertainty_column = resolve_column(
        merged,
        [
            "prediction_uncertainty",
            "ensemble_std",
            "uncertainty",
        ],
        required=False,
    )
    if uncertainty_column:
        merged["uncertainty_component"] = robust_scale(
            merged[uncertainty_column],
            higher_is_better=False,
        )
        uncertainty_available = True
    else:
        merged["uncertainty_component"] = np.nan
        uncertainty_available = False

    merged["pose_quality_component"] = merged.apply(
        pose_quality_score,
        axis=1,
    )
    merged["interaction_component"] = merged.apply(
        interaction_score,
        axis=1,
    )

    pains_column = resolve_column(
        merged,
        ["pains_flag"],
        required=False,
    )
    brenk_column = resolve_column(
        merged,
        ["brenk_flag"],
        required=False,
    )

    if pains_column is None:
        merged["pains_flag"] = False
        pains_column = "pains_flag"

    if brenk_column is None:
        merged["brenk_flag"] = False
        brenk_column = "brenk_flag"

    merged["structural_alert_component"] = merged.apply(
        lambda row: alert_score(
            row[pains_column],
            row[brenk_column],
        ),
        axis=1,
    )

    availability_column = resolve_column(
        merged,
        ["availability"],
        required=False,
    )

    if availability_column:
        merged["supplier_component"] = merged[
            availability_column
        ].apply(availability_score)
    else:
        merged["supplier_component"] = 0.5

    # Use the dedicated, unambiguous structure column for diversity.
    if "consensus_smiles" not in merged.columns:
        raise KeyError(
            "consensus_smiles is missing after input merging."
        )

    merged = calculate_diversity(
        merged,
        smiles_column="consensus_smiles",
    )

    default_score, effective_weights = weighted_score(
        merged,
        DEFAULT_WEIGHTS,
    )
    merged["weighted_consensus_score"] = default_score

    for component in DEFAULT_WEIGHTS:
        merged[f"effective_weight_{component}"] = (
            effective_weights.get(component, 0.0)
        )

    scenario_ranks = pd.DataFrame(
        {"compound_id": merged["compound_id"]}
    )

    scenario_weight_rows = []

    for scenario_name, weights in SCENARIOS.items():
        scenario_score, _ = weighted_score(
            merged,
            weights,
        )

        merged[f"{scenario_name}_score"] = scenario_score

        scenario_ranks[
            f"{scenario_name}_rank"
        ] = (
            scenario_score.rank(
                method="first",
                ascending=False,
            )
            .astype("Int64")
        )

        for component, weight in weights.items():
            scenario_weight_rows.append(
                {
                    "scenario": scenario_name,
                    "component": component,
                    "weight": weight,
                }
            )

    merged["consensus_rank"] = (
        merged["weighted_consensus_score"]
        .rank(method="first", ascending=False)
        .astype(int)
    )

    merged = merged.sort_values(
        [
            "consensus_rank",
            "predicted_pki",
            "vina_score_kcal_mol",
            "compound_id",
        ],
        ascending=[True, False, True, True],
    ).reset_index(drop=True)

    rank_columns = [
        column
        for column in scenario_ranks.columns
        if column.endswith("_rank")
    ]

    scenario_ranks["mean_rank"] = scenario_ranks[
        rank_columns
    ].mean(axis=1)
    scenario_ranks["median_rank"] = scenario_ranks[
        rank_columns
    ].median(axis=1)
    scenario_ranks["rank_std"] = scenario_ranks[
        rank_columns
    ].std(axis=1)
    scenario_ranks["top10_frequency"] = (
        scenario_ranks[rank_columns].le(10).mean(axis=1)
    )
    scenario_ranks["top5_frequency"] = (
        scenario_ranks[rank_columns].le(5).mean(axis=1)
    )

    rank_correlation = scenario_ranks[
        rank_columns
    ].corr(method="spearman")

    complete_mask = (
        merged["predicted_pki"].notna()
        & merged["vina_score_kcal_mol"].notna()
    )

    incomplete = merged.loc[~complete_mask].copy()
    complete = merged.loc[complete_mask].copy()

    in_domain = complete[
        complete.get(
            "applicability_domain_status",
            pd.Series("", index=complete.index),
        ).astype(str).eq("in_domain")
    ]

    high_confidence = in_domain.sort_values(
        [
            "weighted_consensus_score",
            "structural_alert_component",
        ],
        ascending=[False, False],
    )

    top10_diverse_rows = []
    selected_fps = []

    for _, row in complete.sort_values(
        "weighted_consensus_score",
        ascending=False,
    ).iterrows():
        mol = Chem.MolFromSmiles(
            str(row["consensus_smiles"])
        )
        if mol is None:
            continue

        fp = AllChem.GetMorganGenerator(
            radius=2,
            fpSize=2048,
        ).GetFingerprint(mol)

        if selected_fps:
            max_similarity = max(
                DataStructs.BulkTanimotoSimilarity(
                    fp,
                    selected_fps,
                )
            )
            if max_similarity > 0.70:
                continue

        selected_fps.append(fp)
        top10_diverse_rows.append(row)

        if len(top10_diverse_rows) == 10:
            break

    top10_diverse = pd.DataFrame(top10_diverse_rows)

    output_files = {
        "final_ranked_candidates.csv": merged,
        "component_scores.csv": merged[
            [
                "compound_id",
                "predicted_pki_component",
                "docking_component",
                "applicability_component",
                "pose_quality_component",
                "interaction_component",
                "uncertainty_component",
                "structural_alert_component",
                "diversity_component",
                "supplier_component",
                "weighted_consensus_score",
                "consensus_rank",
            ]
        ],
        "top15_review_pool.csv": complete.head(15),
        "top10_consensus_candidates.csv": complete.head(10),
        "top5_high_confidence_candidates.csv": high_confidence.head(5),
        "top10_diverse_candidates.csv": top10_diverse,
        "structural_alert_review.csv": merged[
            merged[pains_column].apply(as_bool)
            | merged[brenk_column].apply(as_bool)
        ],
        "out_of_domain_candidates.csv": merged[
            merged.get(
                "applicability_domain_status",
                pd.Series("", index=merged.index),
            ).astype(str).eq("out_of_domain")
        ],
        "incomplete_evidence_candidates.csv": incomplete,
        "ranking_sensitivity.csv": scenario_ranks,
        "scenario_weights.csv": pd.DataFrame(
            scenario_weight_rows
        ),
    }

    for filename, frame in output_files.items():
        frame.to_csv(
            output_dir / filename,
            index=False,
        )

    rank_correlation.to_csv(
        output_dir / "rank_correlation.csv"
    )

    checksum_map = {
        filename: sha256(output_dir / filename)
        for filename in output_files
    }
    checksum_map["rank_correlation.csv"] = sha256(
        output_dir / "rank_correlation.csv"
    )

    manifest = {
        "schema_version": "1.0",
        "study_id": args.study_id,
        "candidate_count": int(len(merged)),
        "complete_evidence_count": int(
            complete_mask.sum()
        ),
        "incomplete_evidence_count": int(
            (~complete_mask).sum()
        ),
        "in_domain_count": int(
            merged.get(
                "applicability_domain_status",
                pd.Series("", index=merged.index),
            ).astype(str).eq("in_domain").sum()
        ),
        "borderline_count": int(
            merged.get(
                "applicability_domain_status",
                pd.Series("", index=merged.index),
            ).astype(str).eq("borderline").sum()
        ),
        "out_of_domain_count": int(
            merged.get(
                "applicability_domain_status",
                pd.Series("", index=merged.index),
            ).astype(str).eq("out_of_domain").sum()
        ),
        "uncertainty_available": uncertainty_available,
        "default_weights": DEFAULT_WEIGHTS,
        "sensitivity_scenarios": list(
            SCENARIOS
        ),
        "top10_candidate_ids": (
            complete.head(10)["compound_id"].tolist()
        ),
        "top5_high_confidence_ids": (
            high_confidence.head(5)[
                "compound_id"
            ].tolist()
        ),
        "input_paths": {
            "screening": str(args.screening),
            "docking": str(args.docking),
            "pose_summary": str(args.pose_summary),
        },
        "input_checksums": {
            str(path): sha256(path)
            for path in paths
        },
        "output_checksums": checksum_map,
        "git_commit": git_commit(),
        "generated_timestamp": datetime.now(
            timezone.utc
        ).isoformat(),
        "warnings": [
            "Consensus ranking is computational prioritization only.",
            "Vina scores are approximate ranking values.",
            "PAINS and Brenk flags are review cautions, not toxicity conclusions.",
            "Supplier stock must be reconfirmed before purchase.",
        ],
        "phase4_complete": bool(
            len(merged) == 30
            and complete_mask.sum() > 0
            and len(SCENARIOS) >= 5
        ),
        "phase5_ready": bool(
            len(high_confidence) > 0
            and len(merged) == 30
        ),
    }

    manifest_path = output_dir / "consensus_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    schema = {
        "columns": {
            column: str(dtype)
            for column, dtype in merged.dtypes.items()
        }
    }
    (output_dir / "consensus_schema.json").write_text(
        json.dumps(schema, indent=2),
        encoding="utf-8",
    )

    checksums = {
        path.name: sha256(path)
        for path in output_dir.iterdir()
        if path.is_file()
        and path.name != "consensus_checksums.json"
    }
    (
        output_dir / "consensus_checksums.json"
    ).write_text(
        json.dumps(checksums, indent=2),
        encoding="utf-8",
    )

    summary = f"""# DockForgeAI Phase 4 Consensus Ranking

- Study ID: {args.study_id}
- Candidates ranked: {len(merged)}
- Complete evidence: {int(complete_mask.sum())}
- Incomplete evidence: {int((~complete_mask).sum())}
- In-domain candidates: {manifest['in_domain_count']}
- Borderline candidates: {manifest['borderline_count']}
- Out-of-domain candidates: {manifest['out_of_domain_count']}
- Uncertainty available: {uncertainty_available}
- Phase 4 complete: {manifest['phase4_complete']}
- Phase 5 ready for expert review: {manifest['phase5_ready']}

This output is computational prioritization only and does not establish
experimental inhibition, binding, efficacy, or drug-candidate status.
"""
    (output_dir / "consensus_summary.md").write_text(
        summary,
        encoding="utf-8",
    )

    print(json.dumps(manifest, indent=2))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
