"""Generate publication-readable report figures with matplotlib."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Callable

import matplotlib.pyplot as plt

LOGGER = logging.getLogger(__name__)


def _floats(rows: list[dict[str, Any]], key: str) -> list[float]:
    values = []
    for row in rows:
        try:
            values.append(float(row[key]))
        except (KeyError, TypeError, ValueError):
            return []
    return values


def _save(output: Path, draw: Callable[[Any], None]) -> Path | None:
    try:
        output.parent.mkdir(parents=True, exist_ok=True)
        figure, axis = plt.subplots(figsize=(8, 4.8), constrained_layout=True)
        draw(axis)
        figure.savefig(output, dpi=220, bbox_inches="tight")
        return output
    except (OSError, ValueError, TypeError) as exc:
        LOGGER.warning("Skipping figure %s: %s", output.name, exc)
        return None
    finally:
        plt.close("all")


def generate_report_figures(data: dict[str, Any], output_dir: Path) -> list[Path]:
    """Generate every figure supported by the available report sections."""

    outputs: list[Path] = []
    poses = data.get("redocking", {}).get("poses", [])
    ranks = _floats(poses, "pose_rank")
    affinity = _floats(poses, "affinity_kcal_mol")
    rmsd = _floats(poses, "crystal_in_place_rmsd_angstrom")

    if ranks and affinity:
        result = _save(output_dir / "docking_affinity_by_pose.png", lambda ax: (
            ax.bar(ranks, affinity), ax.set_title("Docking affinity by pose"),
            ax.set_xlabel("Pose rank"), ax.set_ylabel("Vina score (kcal/mol)"),
            ax.set_xticks(ranks),
        ))
        if result: outputs.append(result)
    if ranks and rmsd:
        threshold = float(data.get("redocking", {}).get("threshold_angstrom") or 2.0)
        def draw_rmsd(ax: Any) -> None:
            ax.plot(ranks, rmsd, marker="o")
            ax.axhline(threshold, linestyle="--", color="red", label=f"{threshold:g} Å threshold")
            ax.set(title="Crystallographic RMSD by pose", xlabel="Pose rank", ylabel="Heavy-atom RMSD (Å)")
            ax.set_xticks(ranks); ax.legend()
        result = _save(output_dir / "crystallographic_rmsd_by_pose.png", draw_rmsd)
        if result: outputs.append(result)

    interactions = data.get("interactions", {})
    counts = interactions.get("counts_by_type", {})
    if counts:
        result = _save(output_dir / "interaction_counts.png", lambda ax: (
            ax.bar(list(counts), list(counts.values())), ax.set_title("Interaction counts by type"),
            ax.set_xlabel("Interaction type"), ax.set_ylabel("Contact count"), ax.tick_params(axis="x", rotation=25),
        ))
        if result: outputs.append(result)
    residues = interactions.get("top_residues", [])
    if residues:
        result = _save(output_dir / "top_contacted_residues.png", lambda ax: (
            ax.bar([r[0] for r in residues], [r[1] for r in residues]),
            ax.set_title("Top contacted residues"), ax.set_xlabel("Residue"),
            ax.set_ylabel("Contact count"), ax.tick_params(axis="x", rotation=35),
        ))
        if result: outputs.append(result)

    screening = data.get("screening", {})
    screening_rows = screening.get("ranking", [])
    scores = _floats(screening_rows, "best_affinity")
    if scores:
        ids = [str(row.get("compound_id", "?")) for row in screening_rows]
        result = _save(output_dir / "screening_score_ranking.png", lambda ax: (
            ax.bar(ids, scores), ax.set_title("Demonstration screening ranking"),
            ax.set_xlabel("Compound"), ax.set_ylabel("Best Vina score (kcal/mol)"),
        ))
        if result: outputs.append(result)
    workflow = {
        "Accepted": screening.get("accepted_count", 0), "Rejected": screening.get("rejected_count", 0),
        "Succeeded": screening.get("success_count", 0), "Failed": screening.get("failure_count", 0),
    }
    if any(workflow.values()):
        result = _save(output_dir / "screening_workflow_counts.png", lambda ax: (
            ax.bar(list(workflow), list(workflow.values())), ax.set_title("Screening workflow counts"),
            ax.set_xlabel("Outcome"), ax.set_ylabel("Compound count"),
        ))
        if result: outputs.append(result)

    ranking = data.get("prioritization", {}).get("ranking", [])
    consensus = _floats(ranking, "consensus_score")
    if consensus:
        ids = [str(row.get("compound_id", "?")) + f" (#{i + 1})" for i, row in enumerate(ranking)]
        result = _save(output_dir / "consensus_prioritization.png", lambda ax: (
            ax.bar(ids, consensus), ax.set_title("Consensus prioritization ranking"),
            ax.set_xlabel("Compound"), ax.set_ylabel("Consensus score"), ax.tick_params(axis="x", rotation=25),
        ))
        if result: outputs.append(result)
        top = ranking[0]
        descriptor_keys = [key for key in ("qed", "interaction", "ligand_efficiency") if top.get(key) not in (None, "")]
        try:
            descriptor_values = [float(top[key]) for key in descriptor_keys]
        except (TypeError, ValueError):
            descriptor_values = []
        if descriptor_values:
            result = _save(output_dir / "top_hit_descriptor_profile.png", lambda ax: (
                ax.bar(descriptor_keys, descriptor_values), ax.set_title("Top-hit descriptor profile"),
                ax.set_xlabel("Descriptor"), ax.set_ylabel("Reported value"), ax.tick_params(axis="x", rotation=20),
            ))
            if result: outputs.append(result)
    return outputs
