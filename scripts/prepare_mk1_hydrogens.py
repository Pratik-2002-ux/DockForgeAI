"""Prepare crystallographic MK1 with explicit hydrogens."""

from __future__ import annotations

import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


from src.ligands.hydrogenator import (  # noqa: E402
    LigandHydrogenationError,
    prepare_mk1_hydrogens,
)


def main() -> int:
    """Run MK1 hydrogen preparation."""

    try:
        report, report_path = prepare_mk1_hydrogens()
    except (
        FileNotFoundError,
        ValueError,
        KeyError,
        LigandHydrogenationError,
    ) as error:
        print(f"\nMK1 hydrogenation failed: {error}")
        return 1

    print("\nMK1 hydrogen preparation")
    print("-" * 60)
    print(f"Heavy atoms:                     {report['heavy_atoms']}")
    print(f"Total atoms:                     {report['total_atoms']}")
    print(f"Explicit hydrogens:              {report['explicit_hydrogens']}")
    print(f"Formal charge:                   {report['formal_charge']}")
    print(f"Hydrogen optimization:           {report['hydrogen_optimization']}")
    print(f"Heavy-atom RMSD:                 {report['heavy_atom_rmsd_angstrom']} Å")
    print(
        "Maximum heavy-atom displacement: "
        f"{report['maximum_heavy_atom_displacement_angstrom']} Å"
    )
    print(
        "Heavy coordinates preserved:     "
        f"{report['heavy_atom_coordinates_preserved']}"
    )
    print(f"\nPrepared ligand: {report['output_file']}")
    print(f"Report:          {report_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
