"""Create Vina configuration files for screening and refinement."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import yaml

logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SETTINGS_PATH = PROJECT_ROOT / "config" / "settings.yaml"


class ScreeningConfigError(RuntimeError):
    """Raised when screening configuration cannot be generated."""


def load_settings(settings_path: Path = SETTINGS_PATH) -> dict[str, Any]:
    """Load project settings from YAML."""

    with settings_path.open("r", encoding="utf-8") as handle:
        settings = yaml.safe_load(handle)
    if not isinstance(settings, dict):
        raise ValueError("Settings file must contain a YAML dictionary")
    return settings


def write_configs(output_dir: Path) -> dict[str, Path]:
    """Write screening and refinement configs based on the validated redocking box."""

    settings = load_settings()
    output_dir.mkdir(parents=True, exist_ok=True)
    box_path = PROJECT_ROOT / settings["directories"]["docking_boxes"] / "1HSG_MK1_box.json"
    if not box_path.exists():
        raise ScreeningConfigError(f"Docking box not found: {box_path}")

    box = json.loads(box_path.read_text(encoding="utf-8"))
    screening_config = output_dir / "screening_config.txt"
    refinement_config = output_dir / "refinement_config.txt"

    screening_config.write_text(
        "\n".join(
            [
                "receptor = data/processed/pdbqt/receptors/1HSG_chains-A-B_receptor.pdbqt",
                "ligand = {ligand}",
                f"center_x = {box['center']['x']}",
                f"center_y = {box['center']['y']}",
                f"center_z = {box['center']['z']}",
                f"size_x = {box['size']['x']}",
                f"size_y = {box['size']['y']}",
                f"size_z = {box['size']['z']}",
                "exhaustiveness = 4",
                "num_modes = 3",
                "cpu = 4",
                "seed = 2026",
                "energy_range = 3",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    refinement_config.write_text(
        "\n".join(
            [
                "receptor = data/processed/pdbqt/receptors/1HSG_chains-A-B_receptor.pdbqt",
                "ligand = {ligand}",
                f"center_x = {box['center']['x']}",
                f"center_y = {box['center']['y']}",
                f"center_z = {box['center']['z']}",
                f"size_x = {box['size']['x']}",
                f"size_y = {box['size']['y']}",
                f"size_z = {box['size']['z']}",
                "exhaustiveness = 16",
                "num_modes = 10",
                "cpu = 4",
                "seed = 2026",
                "energy_range = 4",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    return {"screening_config": screening_config, "refinement_config": refinement_config}
