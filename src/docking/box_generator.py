"""Generate a docking box for AutoDock Vina from a crystallographic ligand."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import yaml
from rdkit import Chem


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SETTINGS_PATH = PROJECT_ROOT / "config" / "settings.yaml"

DEFAULT_PADDING = 5.0
MIN_BOX_DIMENSION = 20.0
MAX_BOX_DIMENSION = 30.0


class DockingBoxError(RuntimeError):
    """Raised when a docking box cannot be generated."""


def load_settings(settings_path: Path = SETTINGS_PATH) -> dict[str, Any]:
    """Load project settings from YAML."""

    if not settings_path.exists():
        raise FileNotFoundError(f"Settings file not found: {settings_path}")

    with settings_path.open("r", encoding="utf-8") as handle:
        settings = yaml.safe_load(handle)

    if not isinstance(settings, dict):
        raise ValueError("Settings file must contain a YAML dictionary.")

    return settings


def _load_heavy_atom_coordinates(sdf_path: Path) -> list[tuple[float, float, float]]:
    """Load heavy-atom coordinates from an SDF file."""

    if not sdf_path.exists():
        raise FileNotFoundError(f"SDF file not found: {sdf_path}")

    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=True)
    molecule = next((mol for mol in supplier if mol is not None), None)
    if molecule is None:
        raise DockingBoxError(f"Could not read molecule from {sdf_path}")

    conformer = molecule.GetConformer()
    coordinates: list[tuple[float, float, float]] = []
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        position = conformer.GetAtomPosition(atom.GetIdx())
        coordinates.append((float(position.x), float(position.y), float(position.z)))

    if not coordinates:
        raise DockingBoxError("No heavy-atom coordinates found in the ligand SDF")

    return coordinates


def calculate_docking_box(
    ligand_sdf: Path | None = None,
    padding: float = DEFAULT_PADDING,
) -> dict[str, Any]:
    """Calculate a Vina docking box around the crystallographic ligand."""

    settings = load_settings()
    ligand_sdf = ligand_sdf or (
        PROJECT_ROOT / settings["directories"]["crystallographic_ligands"] / "1HSG_MK1_crystallographic_corrected.sdf"
    )

    coordinates = _load_heavy_atom_coordinates(ligand_sdf)
    xs = [coord[0] for coord in coordinates]
    ys = [coord[1] for coord in coordinates]
    zs = [coord[2] for coord in coordinates]

    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    min_z, max_z = min(zs), max(zs)

    center_x = (min_x + max_x) / 2.0
    center_y = (min_y + max_y) / 2.0
    center_z = (min_z + max_z) / 2.0

    span_x = max_x - min_x
    span_y = max_y - min_y
    span_z = max_z - min_z

    size_x = max(MIN_BOX_DIMENSION, span_x + 2.0 * padding)
    size_y = max(MIN_BOX_DIMENSION, span_y + 2.0 * padding)
    size_z = max(MIN_BOX_DIMENSION, span_z + 2.0 * padding)

    warnings: list[str] = []
    for dimension, size in (("x", size_x), ("y", size_y), ("z", size_z)):
        if size > MAX_BOX_DIMENSION:
            warnings.append(f"Box dimension {dimension} is {size:.2f} A, exceeding the recommended 30 A")

    return {
        "ligand_sdf": str(ligand_sdf),
        "padding_angstrom": round(padding, 3),
        "bounds": {
            "min_x": round(min_x, 3),
            "max_x": round(max_x, 3),
            "min_y": round(min_y, 3),
            "max_y": round(max_y, 3),
            "min_z": round(min_z, 3),
            "max_z": round(max_z, 3),
        },
        "center": {
            "x": round(center_x, 3),
            "y": round(center_y, 3),
            "z": round(center_z, 3),
        },
        "size": {
            "x": round(size_x, 3),
            "y": round(size_y, 3),
            "z": round(size_z, 3),
        },
        "span": {
            "x": round(span_x, 3),
            "y": round(span_y, 3),
            "z": round(span_z, 3),
        },
        "warnings": warnings,
    }


def save_box_outputs(box_report: dict[str, Any], settings: dict[str, Any] | None = None) -> tuple[Path, Path]:
    """Save box JSON and Vina config files."""

    settings = settings or load_settings()
    box_directory = PROJECT_ROOT / settings["directories"]["docking_boxes"]
    config_directory = PROJECT_ROOT / settings["directories"]["docking_configs"]
    box_directory.mkdir(parents=True, exist_ok=True)
    config_directory.mkdir(parents=True, exist_ok=True)

    box_json_path = box_directory / "1HSG_MK1_box.json"
    config_path = config_directory / "1HSG_MK1_redocking.txt"

    box_json_path.write_text(json.dumps(box_report, indent=2), encoding="utf-8")

    vina_config = (
        "receptor = data/processed/pdbqt/receptors/1HSG_chains-A-B_receptor.pdbqt\n"
        "ligand = data/processed/pdbqt/ligands/1HSG_MK1.pdbqt\n\n"
        f"center_x = {box_report['center']['x']:.3f}\n"
        f"center_y = {box_report['center']['y']:.3f}\n"
        f"center_z = {box_report['center']['z']:.3f}\n\n"
        f"size_x = {box_report['size']['x']:.3f}\n"
        f"size_y = {box_report['size']['y']:.3f}\n"
        f"size_z = {box_report['size']['z']:.3f}\n\n"
        "exhaustiveness = 8\n"
        "num_modes = 10\n"
        "energy_range = 4\n"
        "cpu = 4\n"
        "seed = 2026\n"
    )
    config_path.write_text(vina_config, encoding="utf-8")

    return box_json_path, config_path


def generate_docking_box() -> dict[str, Any]:
    """Generate the docking box artifacts."""

    settings = load_settings()
    box_report = calculate_docking_box()
    box_json_path, config_path = save_box_outputs(box_report, settings)
    box_report["box_json_path"] = str(box_json_path)
    box_report["vina_config_path"] = str(config_path)
    return box_report


if __name__ == "__main__":
    generate_docking_box()
