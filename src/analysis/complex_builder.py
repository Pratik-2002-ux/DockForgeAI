"""Build protein–ligand complex PDB files from receptor and ligand inputs."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import yaml
from rdkit import Chem

logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SETTINGS_PATH = PROJECT_ROOT / "config" / "settings.yaml"


class ComplexBuilderError(RuntimeError):
    """Raised when a protein–ligand complex cannot be built."""


def load_settings(settings_path: Path = SETTINGS_PATH) -> dict[str, Any]:
    """Load project settings from YAML."""

    with settings_path.open("r", encoding="utf-8") as handle:
        settings = yaml.safe_load(handle)
    if not isinstance(settings, dict):
        raise ValueError("Settings file must contain a YAML dictionary")
    return settings


def build_complex(receptor_path: Path, ligand_path: Path, output_path: Path, label: str) -> dict[str, Any]:
    """Create a simple complex PDB by appending ligand coordinates to the receptor PDB."""

    if not receptor_path.exists():
        raise ComplexBuilderError(f"Receptor PDB not found: {receptor_path}")
    if not ligand_path.exists():
        raise ComplexBuilderError(f"Ligand input not found: {ligand_path}")

    receptor_text = receptor_path.read_text(encoding="utf-8")
    if not receptor_text.strip():
        raise ComplexBuilderError("Receptor PDB is empty")

    receptor_lines = [line for line in receptor_text.splitlines() if line.startswith(("ATOM", "HETATM", "TER"))]
    if not receptor_lines:
        raise ComplexBuilderError("Receptor PDB contains no atom records")

    receptor_chains = sorted({line[21] for line in receptor_lines if len(line) > 21})
    if not receptor_chains:
        raise ComplexBuilderError("Receptor PDB does not preserve chain identifiers")

    ligand_molecule = None
    if ligand_path.suffix.lower() == ".sdf":
        supplier = Chem.SDMolSupplier(str(ligand_path), removeHs=False, sanitize=True)
        molecules = [mol for mol in supplier if mol is not None]
        if not molecules:
            raise ComplexBuilderError("Ligand SDF has no readable molecules")
        ligand_molecule = molecules[0]
        if ligand_molecule.GetNumConformers() == 0:
            raise ComplexBuilderError("Ligand SDF has no conformer")
    else:
        raise ComplexBuilderError("Only SDF ligand input is supported for complex building")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    lines = list(receptor_lines)
    lines.append("TER")

    conformer = ligand_molecule.GetConformer()
    for atom_index in range(ligand_molecule.GetNumAtoms()):
        atom = ligand_molecule.GetAtomWithIdx(atom_index)
        position = conformer.GetAtomPosition(atom_index)
        lines.append(
            f"ATOM  {atom_index + 1:5d} {atom.GetSymbol():<4} LIG L   1    {position.x:8.3f}{position.y:8.3f}{position.z:8.3f}  1.00  0.00           {atom.GetSymbol()}"
        )
    lines.append("END")
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    metadata = {
        "receptor_path": str(receptor_path),
        "ligand_path": str(ligand_path),
        "output_path": str(output_path),
        "label": label,
        "contains_receptor": True,
        "contains_ligand": True,
        "receptor_chains": receptor_chains,
    }
    output_path.with_suffix(output_path.suffix + ".json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return metadata
