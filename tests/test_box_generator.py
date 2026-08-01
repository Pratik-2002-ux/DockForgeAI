import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.docking.box_generator import (
    calculate_docking_box,
    load_settings,
)


def test_box_center_is_within_ligand_bounds() -> None:
    settings = load_settings()
    ligand_path = (
        PROJECT_ROOT
        / settings["directories"]["crystallographic_ligands"]
        / "1HSG_MK1_crystallographic_corrected.sdf"
    )

    report = calculate_docking_box(ligand_path)

    assert report["center"]["x"] >= report["bounds"]["min_x"]
    assert report["center"]["x"] <= report["bounds"]["max_x"]
    assert report["center"]["y"] >= report["bounds"]["min_y"]
    assert report["center"]["y"] <= report["bounds"]["max_y"]
    assert report["center"]["z"] >= report["bounds"]["min_z"]
    assert report["center"]["z"] <= report["bounds"]["max_z"]
    assert report["size"]["x"] >= 20.0
    assert report["size"]["y"] >= 20.0
    assert report["size"]["z"] >= 20.0
