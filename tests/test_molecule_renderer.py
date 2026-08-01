from pathlib import Path
import pytest
from src.visualization.molecule_renderer import MoleculeRenderError, molecule_from_smiles, render_molecule_bytes, save_molecule_images

def test_valid_smiles_renders_png_and_svg() -> None:
    molecule = molecule_from_smiles("CCO")
    assert render_molecule_bytes(molecule, "ethanol", image_format="png").startswith(b"\x89PNG")
    assert b"<svg" in render_molecule_bytes(molecule, "ethanol", image_format="svg")

def test_invalid_smiles_has_clear_error() -> None:
    with pytest.raises(MoleculeRenderError, match="Invalid SMILES"): molecule_from_smiles("not[valid")

def test_save_molecule_images(tmp_path: Path) -> None:
    outputs = save_molecule_images(molecule_from_smiles("CC"), tmp_path / "compound", "compound")
    assert len(outputs) == 2 and all(path.stat().st_size > 0 for path in outputs)
