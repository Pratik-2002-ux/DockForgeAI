"""Validate the RCSB MK1 reference ligand."""

from __future__ import annotations

import sys
from pathlib import Path


# Add DockForgeAI project root to Python path
PROJECT_ROOT = Path(__file__).resolve().parents[1]

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


from src.ligands.validator import (  # noqa: E402
    save_validation_report,
    validate_ligand,
)


def main() -> int:
    """Validate MK1 and print its molecular properties."""

    ligand_path = (
        PROJECT_ROOT
        / "data"
        / "raw"
        / "ligands"
        / "MK1_ideal.sdf"
    )

    if not ligand_path.exists():
        print(f"Error: ligand file not found: {ligand_path}")
        return 1

    try:
        report = validate_ligand(ligand_path)
        report_path = save_validation_report(report)
    except Exception as error:
        print(f"MK1 validation failed: {error}")
        return 1

    properties = report["properties"]

    print("\nLigand validation summary")
    print("-" * 60)
    print(f"Ligand:              {report['ligand_name']}")
    print(f"Sanitization:        {report['sanitization_status']}")
    print(f"3D coordinates:      {report['has_3d_coordinates']}")
    print(f"Conformers:          {report['number_of_conformers']}")
    print(f"Molecular formula:   {properties['molecular_formula']}")
    print(f"Molecular weight:    {properties['molecular_weight']}")
    print(f"Exact molecular wt.: {properties['exact_molecular_weight']}")
    print(f"Heavy atoms:         {properties['heavy_atoms']}")
    print(f"Total atoms:         {properties['total_atoms']}")
    print(f"Formal charge:       {properties['formal_charge']}")
    print(f"H-bond donors:       {properties['hydrogen_bond_donors']}")
    print(f"H-bond acceptors:    {properties['hydrogen_bond_acceptors']}")
    print(f"Rotatable bonds:     {properties['rotatable_bonds']}")
    print(f"Ring count:          {properties['rings']}")
    print(f"Aromatic rings:      {properties['aromatic_rings']}")
    print(f"LogP:                {properties['logp']}")
    print(f"TPSA:                {properties['tpsa']}")
    print(f"Canonical SMILES:    {properties['canonical_smiles']}")
    print(f"\nValidation report: {report_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())