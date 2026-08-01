"""Configuration and rendering helpers for the interactive 3D viewer."""

from __future__ import annotations

import importlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SETTINGS_PATH = PROJECT_ROOT / "config" / "settings.yaml"


@dataclass(frozen=True)
class StructureFile:
    """A molecular structure file and its 3Dmol format."""

    path: Path
    format: str


def _project_path(value: str) -> Path:
    """Resolve a configured path while keeping it inside the project."""

    path = (PROJECT_ROOT / value).resolve()
    try:
        path.relative_to(PROJECT_ROOT.resolve())
    except ValueError as exc:
        raise ValueError(f"Viewer path is outside the project: {value}") from exc
    return path


def load_viewer_settings(settings_path: Path = SETTINGS_PATH) -> dict[str, Any]:
    """Load and validate the viewer section of the project settings."""

    try:
        settings = yaml.safe_load(settings_path.read_text(encoding="utf-8")) or {}
    except (OSError, yaml.YAMLError) as exc:
        raise ValueError(f"Could not load viewer settings: {exc}") from exc

    viewer = settings.get("visualization", {}).get("protein_ligand_viewer")
    if not isinstance(viewer, dict):
        raise ValueError("Missing visualization.protein_ligand_viewer settings")
    return viewer


def structure_files(
    settings: dict[str, Any],
) -> tuple[StructureFile, dict[str, StructureFile]]:
    """Resolve the receptor and configured ligand poses."""

    receptor_config = settings.get("receptor", {})
    receptor = StructureFile(
        _project_path(str(receptor_config.get("path", ""))),
        str(receptor_config.get("format", "pdb")),
    )

    pose_files: dict[str, StructureFile] = {}
    for label, pose in settings.get("poses", {}).items():
        if isinstance(pose, dict):
            pose_files[str(label)] = StructureFile(
                _project_path(str(pose.get("path", ""))),
                str(pose.get("format", "pdbqt")),
            )
    return receptor, pose_files


def resolve_best_rmsd_pose(settings: dict[str, Any]) -> StructureFile:
    """Resolve the best-RMSD docked pose using the generated summary."""

    summary_path = _project_path(str(settings["redocking_summary_path"]))
    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        rank = int(summary["min_rmsd_pose_rank"])
    except (OSError, json.JSONDecodeError, KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"Could not determine the best-RMSD pose: {exc}") from exc

    pose_directory = _project_path(str(settings["docking_pose_directory"]))
    return StructureFile(pose_directory / f"pose_{rank:03d}.pdbqt", "pdbqt")


def available_pose_files(settings: dict[str, Any]) -> dict[str, StructureFile]:
    """Return all three user-facing pose choices in display order."""

    _, configured = structure_files(settings)
    poses = {
        "Crystallographic pose": configured["crystallographic"],
        "Rank-1 redocked pose": configured["rank_1_redocked"],
        "Best-RMSD pose": resolve_best_rmsd_pose(settings),
    }
    return poses


def import_py3dmol() -> Any | None:
    """Import py3Dmol when installed, otherwise signal the optional fallback."""

    try:
        return importlib.import_module("py3Dmol")
    except ImportError:
        return None


def build_view(
    py3dmol: Any,
    receptor: StructureFile,
    ligand: StructureFile,
    catalytic_residues: list[str],
    width: int = 1000,
    height: int = 650,
) -> Any:
    """Create a py3Dmol receptor/ligand view with catalytic-site labels."""

    view = py3dmol.view(width=width, height=height)
    view.addModel(receptor.path.read_text(encoding="utf-8"), receptor.format)
    view.setStyle({"model": 0}, {"cartoon": {"color": "spectrum"}})

    for residue in catalytic_residues:
        chain, residue_number = residue.split(":", maxsplit=1)
        selection = {"model": 0, "chain": chain, "resi": int(residue_number)}
        view.addStyle(
            selection,
            {
                "stick": {"color": "#ef4444", "radius": 0.25},
                "sphere": {"color": "#ef4444", "scale": 0.35},
            },
        )
        view.addLabel(
            f"ASP25 ({chain})",
            {
                "fontColor": "white",
                "backgroundColor": "#991b1b",
                "backgroundOpacity": 0.85,
                "borderRadius": 4,
                "inFront": True,
            },
            selection,
        )

    view.addModel(ligand.path.read_text(encoding="utf-8"), ligand.format)
    view.setStyle(
        {"model": 1},
        {"stick": {"colorscheme": "greenCarbon", "radius": 0.22}},
    )
    view.setBackgroundColor("#0f172a")
    view.zoomTo({"model": 1})
    return view

