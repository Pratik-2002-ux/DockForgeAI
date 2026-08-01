"""Run interaction analysis for protein–ligand complexes."""

from __future__ import annotations

import json
import logging
import math
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any

import yaml

logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SETTINGS_PATH = PROJECT_ROOT / "config" / "settings.yaml"


class InteractionEngineError(RuntimeError):
    """Raised when interaction analysis fails."""


def load_settings(settings_path: Path = SETTINGS_PATH) -> dict[str, Any]:
    """Load project settings from YAML."""

    with settings_path.open("r", encoding="utf-8") as handle:
        settings = yaml.safe_load(handle)
    if not isinstance(settings, dict):
        raise ValueError("Settings file must contain a YAML dictionary")
    return settings


def _parse_pdb_atoms(path: Path) -> list[dict[str, Any]]:
    atoms: list[dict[str, Any]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            parts = line.split()
            if len(parts) >= 9:
                atom_name = parts[2]
                residue_name = parts[3]
                chain = parts[4] if len(parts) > 4 else ""
                residue_number = int(parts[5]) if len(parts) > 5 else 0
                x = float(parts[6])
                y = float(parts[7])
                z = float(parts[8])
            else:
                atom_name = line[12:16].strip()
                residue_name = line[17:20].strip()
                chain = line[21]
                residue_number = int(line[22:26].strip() or 0)
                x = float(line[30:38].strip() or 0.0)
                y = float(line[38:46].strip() or 0.0)
                z = float(line[46:54].strip() or 0.0)
            atoms.append(
                {
                    "atom_name": atom_name,
                    "residue_name": residue_name,
                    "chain": chain,
                    "residue_number": residue_number,
                    "x": x,
                    "y": y,
                    "z": z,
                }
            )
    return atoms


def _fallback_interactions(complex_path: Path, pose_label: str, pose_rank: int, settings: dict[str, Any]) -> list[dict[str, Any]]:
    atoms = _parse_pdb_atoms(complex_path)
    protein_atoms = [atom for atom in atoms if atom["residue_name"] != "LIG"]
    ligand_atoms = [atom for atom in atoms if atom["residue_name"] == "LIG"]
    interactions: list[dict[str, Any]] = []
    hb_distance = settings.get("interactions", {}).get("hydrogen_bond_distance_angstrom", 3.5)
    hb_angle = settings.get("interactions", {}).get("hydrogen_bond_angle_degrees", 120.0)
    hydrophobic_distance = settings.get("interactions", {}).get("hydrophobic_distance_angstrom", 4.0)
    salt_bridge_distance = settings.get("interactions", {}).get("salt_bridge_distance_angstrom", 4.0)
    pi_distance = settings.get("interactions", {}).get("pi_stacking_distance_angstrom", 5.5)

    catalytic_residues = {(item.split(":", 1)[0], int(item.split(":", 1)[1])) for item in settings.get("interactions", {}).get("catalytic_residues", [])}

    for protein_atom in protein_atoms:
        for ligand_atom in ligand_atoms:
            dist = math.dist((protein_atom["x"], protein_atom["y"], protein_atom["z"]), (ligand_atom["x"], ligand_atom["y"], ligand_atom["z"]))
            if dist > 5.5:
                continue
            if protein_atom["atom_name"].upper() in {"OD1", "OD2", "OE1", "OE2", "ND1", "NE2", "NZ", "NH1", "NH2", "OG", "OH"} and ligand_atom["atom_name"].upper() in {"O", "N", "S"} and dist <= hb_distance:
                interactions.append(
                    {
                        "interaction_type": "hydrogen_bond",
                        "protein_chain": protein_atom["chain"],
                        "protein_residue_name": protein_atom["residue_name"],
                        "protein_residue_number": protein_atom["residue_number"],
                        "protein_atom": protein_atom["atom_name"],
                        "ligand_atom": ligand_atom["atom_name"],
                        "distance_angstrom": round(dist, 3),
                        "angle_degrees": hb_angle,
                        "source_engine": "fallback",
                        "pose_label": pose_label,
                        "pose_rank": pose_rank,
                        "notes": "geometric prediction",
                    }
                )
            if protein_atom["atom_name"].upper() in {"CB", "CG", "CD", "CE", "CZ", "CA"} and ligand_atom["atom_name"].upper() in {"C", "CL", "Br", "I"} and dist <= hydrophobic_distance:
                interactions.append(
                    {
                        "interaction_type": "hydrophobic",
                        "protein_chain": protein_atom["chain"],
                        "protein_residue_name": protein_atom["residue_name"],
                        "protein_residue_number": protein_atom["residue_number"],
                        "protein_atom": protein_atom["atom_name"],
                        "ligand_atom": ligand_atom["atom_name"],
                        "distance_angstrom": round(dist, 3),
                        "angle_degrees": None,
                        "source_engine": "fallback",
                        "pose_label": pose_label,
                        "pose_rank": pose_rank,
                        "notes": "geometric prediction",
                    }
                )
            if protein_atom["atom_name"].upper() in {"OD1", "OD2", "OE1", "OE2", "NH1", "NH2", "NZ"} and ligand_atom["atom_name"].upper() in {"O", "N"} and dist <= salt_bridge_distance:
                interactions.append(
                    {
                        "interaction_type": "salt_bridge",
                        "protein_chain": protein_atom["chain"],
                        "protein_residue_name": protein_atom["residue_name"],
                        "protein_residue_number": protein_atom["residue_number"],
                        "protein_atom": protein_atom["atom_name"],
                        "ligand_atom": ligand_atom["atom_name"],
                        "distance_angstrom": round(dist, 3),
                        "angle_degrees": None,
                        "source_engine": "fallback",
                        "pose_label": pose_label,
                        "pose_rank": pose_rank,
                        "notes": "geometric prediction",
                    }
                )
            if protein_atom["residue_name"] in {"PHE", "TYR", "TRP", "HIS"} and ligand_atom["atom_name"].upper() in {"C", "N", "O"} and dist <= pi_distance:
                interactions.append(
                    {
                        "interaction_type": "pi_stacking",
                        "protein_chain": protein_atom["chain"],
                        "protein_residue_name": protein_atom["residue_name"],
                        "protein_residue_number": protein_atom["residue_number"],
                        "protein_atom": protein_atom["atom_name"],
                        "ligand_atom": ligand_atom["atom_name"],
                        "distance_angstrom": round(dist, 3),
                        "angle_degrees": None,
                        "source_engine": "fallback",
                        "pose_label": pose_label,
                        "pose_rank": pose_rank,
                        "notes": "geometric prediction",
                    }
                )
            if (protein_atom["chain"], protein_atom["residue_number"]) in catalytic_residues:
                if dist <= 3.5:
                    interactions.append(
                        {
                            "interaction_type": "catalytic_contact",
                            "protein_chain": protein_atom["chain"],
                            "protein_residue_name": protein_atom["residue_name"],
                            "protein_residue_number": protein_atom["residue_number"],
                            "protein_atom": protein_atom["atom_name"],
                            "ligand_atom": ligand_atom["atom_name"],
                            "distance_angstrom": round(dist, 3),
                            "angle_degrees": None,
                            "source_engine": "fallback",
                            "pose_label": pose_label,
                            "pose_rank": pose_rank,
                            "notes": "catalytic residue contact",
                        }
                    )
    return interactions


def analyze_complex(complex_path: Path, pose_label: str, pose_rank: int, engine: str = "auto") -> dict[str, Any]:
    """Analyze a complex using PLIP if present, otherwise the geometric fallback engine."""

    if not complex_path.exists():
        raise InteractionEngineError(f"Complex PDB not found: {complex_path}")

    settings = load_settings()
    if engine not in {"auto", "plip", "fallback"}:
        raise InteractionEngineError(f"Unsupported interaction engine: {engine}")

    command: list[str] = []
    version = "fallback"
    stdout = ""
    stderr = ""
    runtime_seconds = 0.0
    return_code = 0
    interactions: list[dict[str, Any]] = []

    if engine in {"auto", "plip"}:
        plip = shutil.which("plip")
        if plip:
            command = [plip, str(complex_path)]
            start = time.perf_counter()
            completed = subprocess.run(command, capture_output=True, text=True, check=False, timeout=60)
            runtime_seconds = round(time.perf_counter() - start, 3)
            stdout = completed.stdout
            stderr = completed.stderr
            return_code = completed.returncode
            version = "plip"
            if completed.returncode == 0:
                for line in stdout.splitlines():
                    if line.startswith("INTERACTION"):
                        parts = line.split("\t")
                        if len(parts) >= 3:
                            interactions.append(
                                {
                                    "interaction_type": parts[0],
                                    "protein_chain": parts[1],
                                    "protein_residue_name": parts[2],
                                    "protein_residue_number": 0,
                                    "protein_atom": "",
                                    "ligand_atom": "",
                                    "distance_angstrom": None,
                                    "angle_degrees": None,
                                    "source_engine": "plip",
                                    "pose_label": pose_label,
                                    "pose_rank": pose_rank,
                                    "notes": "plip output",
                                }
                            )
        if not interactions and engine == "plip":
            raise InteractionEngineError("PLIP execution failed")

    if not interactions:
        interactions = _fallback_interactions(complex_path, pose_label, pose_rank, settings)
        version = "fallback"

    return {
        "complex_path": str(complex_path),
        "engine": version,
        "command": " ".join(command),
        "version": version,
        "stdout": stdout,
        "stderr": stderr,
        "runtime_seconds": runtime_seconds,
        "return_code": return_code,
        "interactions": interactions,
    }
