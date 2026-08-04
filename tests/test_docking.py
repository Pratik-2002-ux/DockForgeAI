"""Tests for Phase 4 docking execution and parsing."""
from pathlib import Path

import pytest

from src.docking.phase4_runner import DockingError, parse_scores, run_docking, write_best_pose

ROOT = Path(__file__).resolve().parents[1]


def test_parse_computed_vina_scores(tmp_path: Path) -> None:
    path = tmp_path / "poses.pdbqt"; path.write_text("MODEL 1\nREMARK VINA RESULT: -8.1 0.0 0.0\nENDMDL\n")
    assert parse_scores(path)[0]["docking_score_kcal_mol"] == -8.1


def test_extract_best_pose(tmp_path: Path) -> None:
    source = tmp_path / "poses.pdbqt"; source.write_text("MODEL 1\nREMARK VINA RESULT: -8.1 0 0\nENDMDL\nMODEL 2\nREMARK VINA RESULT: -7 1 2\nENDMDL\n")
    output = tmp_path / "best.pdbqt"; write_best_pose(source, output)
    assert "MODEL 2" not in output.read_text()


def test_optional_engine_stops_gracefully(tmp_path: Path) -> None:
    with pytest.raises(DockingError, match="unavailable|not configured"):
        run_docking(tmp_path, tmp_path / "r", tmp_path / "b", tmp_path / "o", ROOT, engine="quickvina2")
