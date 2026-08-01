#!/usr/bin/env python3
"""Prepare the receptor and generate the docking box for Vina."""

from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.docking.box_generator import generate_docking_box
from src.proteins.receptor_preparer import prepare_receptor


def main() -> int:
    """Run the receptor preparation and docking-box workflow."""

    try:
        receptor_report = prepare_receptor()
        box_report = generate_docking_box()
    except Exception as error:  # pragma: no cover - CLI wrapper
        print(f"Preparation failed: {error}")
        return 1

    receptor_path = Path(receptor_report["receptor_pdbqt"])
    box_json_path = Path(box_report["box_json_path"])
    config_path = Path(box_report["vina_config_path"])

    if not receptor_path.exists() or receptor_path.stat().st_size == 0:
        print("Preparation failed: receptor PDBQT missing or empty")
        return 1
    if not box_json_path.exists() or box_json_path.stat().st_size == 0:
        print("Preparation failed: docking-box JSON missing or empty")
        return 1
    if not config_path.exists() or config_path.stat().st_size == 0:
        print("Preparation failed: Vina config missing or empty")
        return 1

    print(
        "Prepared receptor and docking box: "
        f"{receptor_path.name}, {box_json_path.name}, {config_path.name}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
