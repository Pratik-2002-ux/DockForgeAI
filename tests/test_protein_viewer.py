"""Tests for the optional interactive protein-ligand viewer."""

from __future__ import annotations

from unittest.mock import MagicMock

from app.components.protein_viewer import (
    available_pose_files,
    build_view,
    load_viewer_settings,
    structure_files,
)


def test_configured_viewer_files_exist() -> None:
    settings = load_viewer_settings()
    receptor, _ = structure_files(settings)
    poses = available_pose_files(settings)

    assert receptor.path.is_file()
    assert list(poses) == [
        "Crystallographic pose",
        "Rank-1 redocked pose",
        "Best-RMSD pose",
    ]
    assert all(pose.path.is_file() for pose in poses.values())


def test_build_view_applies_required_representations() -> None:
    settings = load_viewer_settings()
    receptor, _ = structure_files(settings)
    ligand = available_pose_files(settings)["Crystallographic pose"]
    view = MagicMock()
    py3dmol = MagicMock()
    py3dmol.view.return_value = view

    assert build_view(
        py3dmol,
        receptor,
        ligand,
        ["A:25", "B:25"],
    ) is view

    view.setStyle.assert_any_call(
        {"model": 0}, {"cartoon": {"color": "spectrum"}}
    )
    view.setStyle.assert_any_call(
        {"model": 1},
        {"stick": {"colorscheme": "greenCarbon", "radius": 0.22}},
    )
    assert view.addLabel.call_count == 2
    assert view.addStyle.call_count == 2
    view.zoomTo.assert_called_once_with({"model": 1})
