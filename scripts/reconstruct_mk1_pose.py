"""Reconstruct the crystallographic MK1 pose with correct bond orders."""

from __future__ import annotations

import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


from src.ligands.pose_reconstructor import (  # noqa: E402
    PoseReconstructionError,
    reconstruct_mk1_pose,
)


def main() -> int:
    """Run MK1 pose reconstruction."""

    try:
        report, report_path = reconstruct_mk1_pose()
    except (
        FileNotFoundError,
        ValueError,
        KeyError,
        PoseReconstructionError,
    ) as error:
        print(f"\nMK1 pose reconstruction failed: {error}")
        return 1

    print("\nCrystallographic pose reconstruction")
    print("-" * 60)
    print(
        f"Template heavy atoms:      "
        f"{report['template_heavy_atoms']}"
    )
    print(
        f"Pose heavy atoms:          "
        f"{report['crystallographic_heavy_atoms']}"
    )
    print(
        f"Output heavy atoms:        "
        f"{report['output_heavy_atoms']}"
    )
    print(
        f"Output bonds:              "
        f"{report['output_bonds']}"
    )
    print(
        f"Molecular formula:         "
        f"{report['molecular_formula']}"
    )
    print(
        f"Formal charge:             "
        f"{report['formal_charge']}"
    )
    print(
        f"Sanitization:              "
        f"{report['sanitization_status']}"
    )
    print(
        f"Coordinate RMSD:           "
        f"{report['coordinate_rmsd_angstrom']} Å"
    )
    print(
        f"Maximum atom displacement: "
        f"{report['maximum_atom_displacement_angstrom']} Å"
    )
    print(
        f"Coordinates preserved:     "
        f"{report['coordinates_preserved']}"
    )
    print(f"\nCorrected pose: {report['output_file']}")
    print(f"Report:         {report_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())