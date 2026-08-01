"""Generate 3D conformers for standardized ligands."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem import AllChem

logger = logging.getLogger(__name__)


class ConformerGeneratorError(RuntimeError):
    """Raised when a conformer cannot be generated."""


def generate_conformers(records: list[dict[str, Any]], output_dir: Path) -> list[dict[str, Any]]:
    """Generate one 3D conformer per ligand and write an SDF file per compound."""

    output_dir.mkdir(parents=True, exist_ok=True)
    generated: list[dict[str, Any]] = []
    for record in records:
        mol = Chem.MolFromSmiles(record["smiles"])
        if mol is None:
            raise ConformerGeneratorError(f"SMILES could not be parsed for {record['compound_id']}")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=0)
        try:
            AllChem.MMFFOptimizeMolecule(mol)
            forcefield = "MMFF"
        except Exception:
            AllChem.UFFOptimizeMolecule(mol)
            forcefield = "UFF"
        output_path = output_dir / f"{record['compound_id']}.sdf"
        writer = Chem.SDWriter(str(output_path))
        writer.write(mol)
        writer.close()
        generated.append({"compound_id": record["compound_id"], "output_path": str(output_path), "forcefield": forcefield})
    return generated
