"""Deterministic RDKit 2D molecule rendering for reports and Streamlit."""

from __future__ import annotations

from io import BytesIO
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem import AllChem, Draw


class MoleculeRenderError(ValueError):
    """Raised when a molecule cannot be parsed or rendered."""


def molecule_from_smiles(smiles: str) -> Chem.Mol:
    """Parse SMILES and generate deterministic 2D coordinates."""

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise MoleculeRenderError(f"Invalid SMILES: {smiles}")
    AllChem.Compute2DCoords(molecule, canonOrient=True)
    return molecule


def molecule_from_sdf(path: Path) -> Chem.Mol:
    """Read the first valid molecule from an SDF file."""

    if not path.is_file():
        raise MoleculeRenderError(f"SDF file not found: {path.name}")
    molecule = next((mol for mol in Chem.SDMolSupplier(str(path), removeHs=False) if mol), None)
    if molecule is None:
        raise MoleculeRenderError(f"No valid molecule in SDF: {path.name}")
    AllChem.Compute2DCoords(molecule, canonOrient=True)
    return molecule


def render_molecule_bytes(
    molecule: Chem.Mol,
    compound_id: str,
    score_label: str | None = None,
    image_format: str = "png",
    size: tuple[int, int] = (500, 400),
) -> bytes:
    """Render a molecule as PNG or SVG bytes."""

    if molecule is None:
        raise MoleculeRenderError("A valid RDKit molecule is required")
    legend = compound_id if not score_label else f"{compound_id} | {score_label}"
    fmt = image_format.lower()
    try:
        if fmt == "png":
            image = Draw.MolToImage(molecule, size=size, legend=legend)
            buffer = BytesIO()
            image.save(buffer, format="PNG")
            return buffer.getvalue()
        if fmt == "svg":
            drawer = Draw.MolDraw2DSVG(*size)
            drawer.DrawMolecule(molecule, legend=legend)
            drawer.FinishDrawing()
            return drawer.GetDrawingText().encode("utf-8")
    except Exception as exc:  # RDKit exposes multiple renderer exception types
        raise MoleculeRenderError(f"Could not render {compound_id}: {exc}") from exc
    raise MoleculeRenderError(f"Unsupported image format: {image_format}")


def save_molecule_images(
    molecule: Chem.Mol,
    output_stem: Path,
    compound_id: str,
    score_label: str | None = None,
) -> list[Path]:
    """Save PNG and SVG renderings and return their paths."""

    output_stem.parent.mkdir(parents=True, exist_ok=True)
    outputs = []
    for fmt in ("png", "svg"):
        output = output_stem.with_suffix(f".{fmt}")
        output.write_bytes(render_molecule_bytes(molecule, compound_id, score_label, fmt))
        outputs.append(output)
    return outputs
