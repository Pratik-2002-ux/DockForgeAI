from pathlib import Path

import pytest

from src.docking.pose_exporter import PoseExportError, export_poses


def test_export_poses_uses_fallback_for_missing_meeko(tmp_path: Path) -> None:
    reference_sdf = tmp_path / "reference.sdf"
    output_sdf = tmp_path / "poses.sdf"
    report_path = tmp_path / "report.json"

    reference_mol = pytest.importorskip("rdkit.Chem").MolFromMolBlock(
        """
  -INDIGO-07202408162D

  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
"""
    )
    if reference_mol is None:
        pytest.skip("RDKit could not build a test molecule")

    from rdkit import Chem

    reference_mol = Chem.AddHs(reference_mol)
    writer = Chem.SDWriter(str(reference_sdf))
    writer.write(reference_mol)
    writer.close()

    pose_pdbqt_1 = tmp_path / "pose1.pdbqt"
    pose_pdbqt_1.write_text(
        """MODEL 1
REMARK VINA RESULT: -10.000 0.000 0.000
ATOM      1  C1  LIG A   1      0.000 0.000 0.000
ATOM      2  C2  LIG A   1      1.000 0.000 0.000
ATOM      3  C3  LIG A   1      0.000 1.000 0.000
ENDMDL
""",
        encoding="utf-8",
    )
    pose_pdbqt_2 = tmp_path / "pose2.pdbqt"
    pose_pdbqt_2.write_text(
        """MODEL 2
REMARK VINA RESULT: -9.000 0.000 0.000
ATOM      1  C1  LIG A   1      0.000 0.000 0.000
ATOM      2  C2  LIG A   1      1.000 0.000 0.000
ATOM      3  C3  LIG A   1      0.000 1.000 0.000
ENDMDL
""",
        encoding="utf-8",
    )

    report = export_poses(
        pose_pdbqts=[pose_pdbqt_1, pose_pdbqt_2],
        output_sdf=output_sdf,
        reference_sdf=reference_sdf,
        report_path=report_path,
    )

    assert report["pose_count"] == 2
    assert report["export_method"] in {"meeko", "fallback_rdkit"}
    assert output_sdf.exists()


def test_export_poses_rejects_mismatched_pose_count(tmp_path: Path) -> None:
    reference_sdf = tmp_path / "reference.sdf"
    output_sdf = tmp_path / "poses.sdf"
    report_path = tmp_path / "report.json"

    from rdkit import Chem

    reference_mol = Chem.MolFromSmiles("CCO")
    if reference_mol is None:
        pytest.skip("RDKit could not build a test molecule")

    writer = Chem.SDWriter(str(reference_sdf))
    writer.write(reference_mol)
    writer.close()

    pose_path = tmp_path / "pose.pdbqt"
    pose_path.write_text(
        """MODEL 1
REMARK VINA RESULT: -10.000 0.000 0.000
ATOM      1  C1  LIG A   1      0.000 0.000 0.000
ATOM      2  C2  LIG A   1      1.000 0.000 0.000
ATOM      3  C3  LIG A   1      0.000 1.000 0.000
ENDMDL
""",
        encoding="utf-8",
    )

    with pytest.raises(PoseExportError):
        export_poses(
            pose_pdbqts=[pose_path],
            output_sdf=output_sdf,
            reference_sdf=reference_sdf,
            report_path=report_path,
        )
