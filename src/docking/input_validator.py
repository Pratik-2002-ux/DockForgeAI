"""Validate docking inputs before running AutoDock Vina."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import yaml
from rdkit import Chem


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SETTINGS_PATH = PROJECT_ROOT / "config" / "settings.yaml"


class DockingInputError(RuntimeError):
    """Raised when a docking input is invalid."""


def load_settings(settings_path: Path = SETTINGS_PATH) -> dict[str, Any]:
    """Load project settings from YAML."""

    if not settings_path.exists():
        raise FileNotFoundError(f"Settings file not found: {settings_path}")

    with settings_path.open("r", encoding="utf-8") as handle:
        settings = yaml.safe_load(handle)

    if not isinstance(settings, dict):
        raise ValueError("Settings file must contain a YAML dictionary.")

    return settings


def resolve_path(settings: dict[str, Any], key: str) -> Path:
    """Resolve a project-relative path from the settings dictionary."""

    value = settings["directories"].get(key)
    if not isinstance(value, str) or not value.strip():
        raise DockingInputError(f"Directory setting '{key}' is not configured")
    return PROJECT_ROOT / value


def _ensure_non_empty(path: Path, description: str) -> None:
    """Verify a file exists and is non-empty."""

    if not path.exists():
        raise DockingInputError(f"{description} not found: {path}")
    if path.stat().st_size == 0:
        raise DockingInputError(f"{description} is empty: {path}")


def _parse_vina_config(config_path: Path) -> dict[str, str]:
    """Parse a simple Vina-style config file."""

    parsed: dict[str, str] = {}
    for line in config_path.read_text(encoding="utf-8").splitlines():
        stripped = line.split("#", 1)[0].strip()
        if not stripped:
            continue
        if "=" not in stripped:
            continue
        key, value = stripped.split("=", 1)
        parsed[key.strip().lower()] = value.strip()
    return parsed


def validate_docking_inputs() -> dict[str, Any]:
    """Validate receptor, ligand, reference, box and config inputs."""

    settings = load_settings()
    receptor_pdbqt = PROJECT_ROOT / settings["directories"]["receptor_pdbqt"] / "1HSG_chains-A-B_receptor.pdbqt"
    ligand_pdbqt = PROJECT_ROOT / settings["directories"]["ligand_pdbqt"] / "1HSG_MK1.pdbqt"
    reference_sdf = PROJECT_ROOT / settings["directories"]["crystallographic_ligands"] / "1HSG_MK1_crystallographic_corrected.sdf"
    docking_config = PROJECT_ROOT / settings["directories"]["docking_configs"] / "1HSG_MK1_redocking.txt"
    box_report_path = PROJECT_ROOT / settings["directories"]["docking_boxes"] / "1HSG_MK1_box.json"
    preparation_report = PROJECT_ROOT / settings["directories"]["preparation_reports"] / "1HSG_receptor_preparation.json"

    for path, description in [
        (receptor_pdbqt, "Receptor PDBQT"),
        (ligand_pdbqt, "Ligand PDBQT"),
        (reference_sdf, "Reference SDF"),
        (docking_config, "Docking config"),
        (box_report_path, "Docking box report"),
        (preparation_report, "Receptor preparation report"),
    ]:
        _ensure_non_empty(path, description)

    receptor_text = receptor_pdbqt.read_text(encoding="utf-8")
    if not any(line.startswith(("ATOM", "HETATM")) for line in receptor_text.splitlines()):
        raise DockingInputError("Receptor PDBQT does not contain atom records")
    if any(token in receptor_text for token in ("ROOT", "BRANCH", "TORSDOF")):
        raise DockingInputError("Receptor PDBQT contains ligand torsion-tree records")

    ligand_text = ligand_pdbqt.read_text(encoding="utf-8")
    if not any(line.startswith(("ATOM", "HETATM")) for line in ligand_text.splitlines()):
        raise DockingInputError("Ligand PDBQT does not contain coordinate records")
    for token in ("ROOT", "ENDROOT", "TORSDOF"):
        if token not in ligand_text:
            raise DockingInputError(f"Ligand PDBQT is missing '{token}'")

    config_values = _parse_vina_config(docking_config)
    required_keys = [
        "receptor",
        "ligand",
        "center_x",
        "center_y",
        "center_z",
        "size_x",
        "size_y",
        "size_z",
        "exhaustiveness",
        "num_modes",
        "energy_range",
        "cpu",
        "seed",
    ]
    missing_keys = [key for key in required_keys if key not in config_values]
    if missing_keys:
        raise DockingInputError(f"Docking config is missing keys: {', '.join(missing_keys)}")

    cpu_value = int(config_values["cpu"])
    if cpu_value > 4:
        raise DockingInputError(f"Configured CPU count exceeds the allowed maximum: {cpu_value}")

    box_report = json.loads(box_report_path.read_text(encoding="utf-8"))
    size = box_report.get("size", {})
    for dimension in ("x", "y", "z"):
        value = float(size[dimension])
        if value < 20.0 or value > 30.0:
            raise DockingInputError(
                f"Docking box dimension {dimension} is outside the allowed range: {value:.3f}"
            )

    supplier = Chem.SDMolSupplier(str(reference_sdf), removeHs=False, sanitize=True)
    molecules = [mol for mol in supplier if mol is not None]
    if len(molecules) != 1:
        raise DockingInputError("Reference SDF should contain a single readable molecule")
    reference_molecule = molecules[0]
    if reference_molecule.GetNumConformers() != 1:
        raise DockingInputError("Reference SDF must contain exactly one conformer")
    heavy_atoms = reference_molecule.GetNumHeavyAtoms()
    if heavy_atoms != 45:
        raise DockingInputError(f"Reference molecule should have 45 heavy atoms, found {heavy_atoms}")

    return {
        "settings": settings,
        "receptor_pdbqt": receptor_pdbqt,
        "ligand_pdbqt": ligand_pdbqt,
        "reference_sdf": reference_sdf,
        "docking_config": docking_config,
        "box_report": box_report,
        "preparation_report": preparation_report,
        "config_values": config_values,
        "reference_molecule": reference_molecule,
    }
