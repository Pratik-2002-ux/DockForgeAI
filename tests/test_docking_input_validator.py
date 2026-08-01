from pathlib import Path

import pytest

from src.docking.input_validator import DockingInputError, validate_docking_inputs


@pytest.fixture()
def project_root() -> Path:
    return Path(__file__).resolve().parents[1]


def test_validate_docking_inputs_accepts_existing_inputs(project_root: Path) -> None:
    validation = validate_docking_inputs()

    assert validation["receptor_pdbqt"].exists()
    assert validation["ligand_pdbqt"].exists()
    assert validation["reference_sdf"].exists()
    assert validation["config_values"]["receptor"].endswith("1HSG_chains-A-B_receptor.pdbqt")
    assert validation["config_values"]["ligand"].endswith("1HSG_MK1.pdbqt")
    assert validation["config_values"]["cpu"] == "4"


def test_validate_docking_inputs_rejects_missing_input(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    import src.docking.input_validator as input_validator

    monkeypatch.setattr(input_validator, "PROJECT_ROOT", tmp_path)

    def fake_load_settings(settings_path: Path = input_validator.SETTINGS_PATH) -> dict[str, object]:
        return {
            "directories": {
                "receptor_pdbqt": "data/processed/pdbqt/receptors",
                "ligand_pdbqt": "data/processed/pdbqt/ligands",
                "crystallographic_ligands": "data/processed/crystallographic_ligands",
                "docking_configs": "config/docking",
                "docking_boxes": "data/processed/docking_boxes",
                "preparation_reports": "data/processed/preparation_reports",
            }
        }

    monkeypatch.setattr(input_validator, "load_settings", fake_load_settings)

    (tmp_path / "data/processed/pdbqt/receptors").mkdir(parents=True, exist_ok=True)
    (tmp_path / "data/processed/pdbqt/ligands").mkdir(parents=True, exist_ok=True)
    (tmp_path / "data/processed/crystallographic_ligands").mkdir(parents=True, exist_ok=True)
    (tmp_path / "config/docking").mkdir(parents=True, exist_ok=True)
    (tmp_path / "data/processed/docking_boxes").mkdir(parents=True, exist_ok=True)
    (tmp_path / "data/processed/preparation_reports").mkdir(parents=True, exist_ok=True)

    (tmp_path / "data/processed/pdbqt/ligands" / "1HSG_MK1.pdbqt").write_text("ATOM 1", encoding="utf-8")
    (tmp_path / "data/processed/crystallographic_ligands" / "1HSG_MK1_crystallographic_corrected.sdf").write_text("", encoding="utf-8")
    (tmp_path / "config/docking" / "1HSG_MK1_redocking.txt").write_text("receptor = x\n", encoding="utf-8")
    (tmp_path / "data/processed/docking_boxes" / "1HSG_MK1_box.json").write_text("{\"size\": {\"x\": 25, \"y\": 25, \"z\": 25}}", encoding="utf-8")
    (tmp_path / "data/processed/preparation_reports" / "1HSG_receptor_preparation.json").write_text("{}", encoding="utf-8")

    with pytest.raises(DockingInputError):
        validate_docking_inputs()
