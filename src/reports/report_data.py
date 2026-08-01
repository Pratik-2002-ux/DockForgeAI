"""Aggregate validated DockForgeAI outputs into report-ready data."""

from __future__ import annotations

import csv
import importlib.metadata
import json
import platform
import subprocess
from collections import Counter
from pathlib import Path
from typing import Any

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SETTINGS_PATH = PROJECT_ROOT / "config" / "settings.yaml"


def _json(path: Path, warnings: list[str], optional: bool = True) -> dict[str, Any]:
    if not path.is_file():
        if optional:
            warnings.append(f"Optional input unavailable: {_relative(path)}")
        return {}
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {}
    except (OSError, json.JSONDecodeError) as exc:
        warnings.append(f"Could not read {_relative(path)}: {exc}")
        return {}


def _csv(path: Path, warnings: list[str]) -> list[dict[str, str]]:
    if not path.is_file():
        warnings.append(f"Optional input unavailable: {_relative(path)}")
        return []
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle))
    except (OSError, csv.Error) as exc:
        warnings.append(f"Could not read {_relative(path)}: {exc}")
        return []


def _relative(path: Path) -> str:
    try:
        return path.resolve().relative_to(PROJECT_ROOT.resolve()).as_posix()
    except ValueError:
        return path.name


def _number(value: Any, default: float | int | None = None) -> float | int | None:
    if isinstance(value, bool):
        return default
    try:
        result = float(value)
    except (TypeError, ValueError):
        return default
    return int(result) if result.is_integer() else result


def _version(distribution: str) -> str:
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return "unavailable"


def _command_version(command: list[str]) -> str:
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=5, check=False)
        return (result.stdout or result.stderr).strip().splitlines()[0]
    except (OSError, subprocess.SubprocessError, IndexError):
        return "unavailable"


def _git_hash(root: Path) -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"], cwd=root,
            capture_output=True, text=True, timeout=5, check=True,
        ).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return "unavailable"


def _sanitize_paths(value: Any, root: Path) -> Any:
    """Replace absolute in-project paths with portable project-relative paths."""

    if isinstance(value, dict):
        return {key: _sanitize_paths(item, root) for key, item in value.items()}
    if isinstance(value, list):
        return [_sanitize_paths(item, root) for item in value]
    if isinstance(value, str) and Path(value).is_absolute():
        try:
            return Path(value).resolve().relative_to(root.resolve()).as_posix()
        except ValueError:
            return Path(value).name
    return value


def collect_report_data(
    project_root: Path = PROJECT_ROOT,
    settings_path: Path | None = None,
    target_id: str = "1HSG",
    ligand_id: str = "MK1",
) -> dict[str, Any]:
    """Collect final outputs safely; missing optional sections become warnings."""

    root = project_root.resolve()
    config_path = settings_path or root / "config" / "settings.yaml"
    try:
        settings = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    except (OSError, yaml.YAMLError) as exc:
        raise ValueError(f"Unable to load settings: {exc}") from exc
    directories = settings.get("directories", {})
    warnings: list[str] = []

    def path(directory: str, name: str) -> Path:
        return root / directories.get(directory, f"data/processed/{directory}") / name

    inspection = _json(path("protein_reports", f"{target_id}_inspection.json"), warnings)
    cleaning = _json(path("preparation_reports", f"{target_id}_cleaning_report.json"), warnings)
    preparation = _json(path("preparation_reports", f"{target_id}_receptor_preparation.json"), warnings)
    redocking = _json(path("docking_results", f"{target_id}_{ligand_id}_redocking_summary.json"), warnings)
    docking_box = _json(path("docking_boxes", f"{target_id}_{ligand_id}_box.json"), warnings)
    interaction_rows = _csv(path("interaction_reports", "interactions.csv"), warnings)
    screening_summary = _json(path("screening_reports", "screening_report.json"), warnings)
    screening_rows = _csv(path("screening_tables", "screening_scores.csv"), warnings)
    prioritization_summary = _json(path("prioritization_reports", "hit_prioritization_report.json"), warnings)
    prioritization_rows = _csv(path("prioritization_reports", "consensus_ranking.csv"), warnings)
    diverse_rows = _csv(path("prioritization_reports", "diverse_top_hits.csv"), warnings)
    pose_rows = _csv(path("docking_tables", f"{target_id}_{ligand_id}_redocking_analysis.csv"), warnings)

    interaction_types = Counter(row.get("interaction_type", "unknown") for row in interaction_rows)
    residues = Counter(
        f"{row.get('protein_residue_name', '')}{row.get('protein_residue_number', '')} ({row.get('protein_chain', '')})"
        for row in interaction_rows
    )
    catalytic = [r for r in interaction_rows if r.get("protein_residue_name") == "ASP" and r.get("protein_residue_number") == "25"]
    pose_counts = Counter(row.get("pose_label", "unknown") for row in interaction_rows)
    successful = [r for r in screening_rows if r.get("status") == "success"]
    best = min(successful, key=lambda r: float(r["best_affinity"])) if successful else {}
    runtime = sum(float(r.get("runtime_seconds") or 0) for r in screening_rows)
    metadata = inspection.get("experimental_metadata", {})
    sources = [
        path("protein_reports", f"{target_id}_inspection.json"),
        path("docking_results", f"{target_id}_{ligand_id}_redocking_summary.json"),
        path("interaction_reports", "interactions.csv"),
        path("screening_tables", "screening_scores.csv"),
        path("prioritization_reports", "consensus_ranking.csv"),
    ]

    return {
        "project": {"name": settings.get("project", {}).get("name", "DockForgeAI"), "version": settings.get("project", {}).get("version", "unknown"), "target_id": target_id, "ligand_id": ligand_id},
        "software": {
            "python": platform.python_version(), "rdkit": _version("rdkit"),
            "vina": _command_version(["vina", "--version"]), "meeko": _version("meeko"),
            "gemmi": _version("gemmi"), "streamlit": _version("streamlit"),
            "py3Dmol": _version("py3Dmol"), "git_commit": _git_hash(root),
        },
        "protein": {
            "pdb_id": inspection.get("structure_id", target_id),
            "title": metadata.get("title", "Unavailable"),
            "experimental_method": ", ".join(metadata.get("experimental_methods", [])) or "Unavailable",
            "resolution_angstrom": (metadata.get("resolution_angstrom") or [None])[0],
            "chains": cleaning.get("selected_chains", [c.get("chain_id") for c in inspection.get("chains", [])]),
            "residue_count": cleaning.get("protein_residues_retained", inspection.get("protein_residues")),
            "removed_waters": cleaning.get("waters_removed", inspection.get("water_molecules")),
            "extracted_ligands": _sanitize_paths(cleaning.get("extracted_ligands", inspection.get("ligands", [])), root),
            "protonation_assumption": preparation.get("protonation_assumption", "Unavailable"),
        },
        "redocking": {
            "affinity_kcal_mol": _number(redocking.get("best_affinity_kcal_mol")),
            "rmsd_angstrom": _number(redocking.get("min_rmsd_angstrom")),
            "pose_rank": _number(redocking.get("min_rmsd_pose_rank")),
            "passing_pose_count": _number(redocking.get("passing_pose_count"), 0),
            "interpretation": redocking.get("interpretation", "Unavailable"),
            "threshold_angstrom": _number(redocking.get("threshold_angstrom"), 2.0),
            "box": _sanitize_paths(docking_box, root), "poses": _sanitize_paths(pose_rows, root),
        },
        "interactions": {
            "total_contacts": len(interaction_rows), "counts_by_type": dict(interaction_types),
            "top_residues": residues.most_common(10), "catalytic_asp25_contacts": len(catalytic),
            "pose_comparison": dict(pose_counts), "rows": interaction_rows,
        },
        "screening": {
            "library_size": screening_summary.get("library_size", len(screening_rows)),
            "accepted_count": screening_summary.get("accepted_count", len(screening_rows)),
            "rejected_count": screening_summary.get("rejected_count", 0),
            "success_count": screening_summary.get("screening_result_count", len(successful)),
            "failure_count": screening_summary.get("screening_failure_count", len(screening_rows) - len(successful)),
            "ranking": screening_rows, "best_demonstration_compound": best,
            "runtime_seconds": round(runtime, 3),
        },
        "prioritization": {
            "ranking": prioritization_rows or prioritization_summary.get("rows", []),
            "diverse_shortlist": diverse_rows,
            "weights": settings.get("prioritization", {}),
            "external_admet_provenance": "Not supplied; RDKit descriptors are not full ADMET predictions.",
        },
        "limitations": [
            "Docking scores are approximate ranking values, not experimental binding free energies.",
            "The three-compound screen demonstrates software operation and does not identify validated inhibitors.",
            "RDKit descriptors are not full ADMET predictions.",
            "All computational predictions require experimental validation.",
        ],
        "files": {"source_inputs": [_relative(p) for p in sources if p.is_file()]},
        "warnings": warnings,
    }
