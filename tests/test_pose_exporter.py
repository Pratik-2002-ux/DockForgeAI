"""Tests for chemically valid export of distinct docked PDBQT poses."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from src.docking.pose_exporter import (
    PoseExportError,
    export_poses,
)


def _write_reference_sdf(path: Path) -> None:
    """Write a deterministic 3D ethanol molecule as an SDF template."""

    molecule = Chem.MolFromSmiles("CCO")

    if molecule is None:
        raise RuntimeError("Could not create ethanol test molecule.")

    molecule = Chem.AddHs(molecule)

    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = 2026

    result = AllChem.EmbedMolecule(
        molecule,
        parameters,
    )

    if result != 0:
        raise RuntimeError("Could not generate ethanol coordinates.")

    AllChem.UFFOptimizeMolecule(molecule)

    molecule = Chem.RemoveHs(molecule)
    molecule.SetProp("_Name", "ethanol_reference")

    writer = Chem.SDWriter(str(path))

    try:
        writer.write(molecule)
    finally:
        writer.close()


def _atom_line(
    serial: int,
    atom_name: str,
    x: float,
    y: float,
    z: float,
    atom_type: str,
) -> str:
    """Create a minimally valid fixed-width PDBQT atom record."""

    return (
        f"ATOM  {serial:5d} {atom_name:>4s} UNL     1"
        f"    {x:8.3f}{y:8.3f}{z:8.3f}"
        f"  1.00  0.00     0.000 {atom_type:<2s}"
    )


def _write_pose_pdbqt(
    path: Path,
    coordinates: list[tuple[float, float, float]],
) -> None:
    """Write a small Meeko-style ethanol PDBQT pose."""

    if len(coordinates) != 3:
        raise ValueError("Ethanol fixture requires three heavy atoms.")

    lines = [
        "REMARK SMILES CCO",
        "REMARK SMILES IDX 1 1 2 2 3 3",
        "ROOT",
        _atom_line(
            1,
            "C1",
            *coordinates[0],
            "C",
        ),
        _atom_line(
            2,
            "C2",
            *coordinates[1],
            "C",
        ),
        _atom_line(
            3,
            "O1",
            *coordinates[2],
            "OA",
        ),
        "ENDROOT",
        "TORSDOF 0",
    ]

    path.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def test_export_poses_preserves_distinct_coordinates(
    tmp_path: Path,
) -> None:
    """Different PDBQT poses must remain different in exported SDF."""

    reference_sdf = tmp_path / "reference.sdf"
    pose_1 = tmp_path / "pose_001.pdbqt"
    pose_2 = tmp_path / "pose_002.pdbqt"
    output_sdf = tmp_path / "poses.sdf"
    report_path = tmp_path / "export_report.json"

    _write_reference_sdf(reference_sdf)

    _write_pose_pdbqt(
        pose_1,
        [
            (0.0, 0.0, 0.0),
            (1.5, 0.0, 0.0),
            (2.8, 0.2, 0.0),
        ],
    )

    _write_pose_pdbqt(
        pose_2,
        [
            (1.0, 1.0, 1.0),
            (2.5, 1.0, 1.0),
            (3.8, 1.2, 1.0),
        ],
    )

    report = export_poses(
        pose_pdbqts=[pose_1, pose_2],
        output_sdf=output_sdf,
        reference_sdf=reference_sdf,
        report_path=report_path,
    )

    assert output_sdf.exists()
    assert output_sdf.stat().st_size > 0
    assert report_path.exists()

    assert report["export_method"] == "rdkit_meeko_smiles_idx"
    assert report["pose_count"] == 2
    assert report["unique_coordinate_sets"] == 2
    assert report["all_poses_identical"] is False

    saved_report = json.loads(
        report_path.read_text(encoding="utf-8")
    )

    assert saved_report["pose_count"] == 2
    assert saved_report["unique_coordinate_sets"] == 2

    supplier = Chem.SDMolSupplier(
        str(output_sdf),
        removeHs=True,
        sanitize=True,
    )

    poses = [
        molecule
        for molecule in supplier
        if molecule is not None
    ]

    assert len(poses) == 2
    assert poses[0].GetNumHeavyAtoms() == 3
    assert poses[1].GetNumHeavyAtoms() == 3

    first_conformer = poses[0].GetConformer()
    second_conformer = poses[1].GetConformer()

    coordinates_differ = any(
        (
            abs(
                first_conformer.GetAtomPosition(index).x
                - second_conformer.GetAtomPosition(index).x
            )
            > 1e-6
            or abs(
                first_conformer.GetAtomPosition(index).y
                - second_conformer.GetAtomPosition(index).y
            )
            > 1e-6
            or abs(
                first_conformer.GetAtomPosition(index).z
                - second_conformer.GetAtomPosition(index).z
            )
            > 1e-6
        )
        for index in range(poses[0].GetNumAtoms())
    )

    assert coordinates_differ


def test_export_poses_rejects_missing_meeko_smiles(
    tmp_path: Path,
) -> None:
    """A PDBQT lacking Meeko chemical mapping must fail clearly."""

    reference_sdf = tmp_path / "reference.sdf"
    invalid_pose = tmp_path / "invalid_pose.pdbqt"
    output_sdf = tmp_path / "poses.sdf"
    report_path = tmp_path / "report.json"

    _write_reference_sdf(reference_sdf)

    invalid_pose.write_text(
        "\n".join(
            [
                "ROOT",
                _atom_line(1, "C1", 0.0, 0.0, 0.0, "C"),
                "ENDROOT",
                "TORSDOF 0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(
        PoseExportError,
        match="REMARK SMILES",
    ):
        export_poses(
            pose_pdbqts=[invalid_pose],
            output_sdf=output_sdf,
            reference_sdf=reference_sdf,
            report_path=report_path,
        )
