import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.proteins.receptor_preparer import (
    load_settings,
    validate_clean_receptor,
)


def test_receptor_contains_required_chains_and_asp25() -> None:
    settings = load_settings()
    receptor_path = (
        PROJECT_ROOT
        / settings["directories"]["prepared_proteins"]
        / "1HSG_chains-A-B_clean.pdb"
    )

    report = validate_clean_receptor(receptor_path)

    assert {"A", "B"}.issubset(report["chains"])
    assert report["catalytic_asp25"]["A"] == 25
    assert report["catalytic_asp25"]["B"] == 25
    assert report["contains_waters"] is False
    assert report["contains_mk1"] is False
