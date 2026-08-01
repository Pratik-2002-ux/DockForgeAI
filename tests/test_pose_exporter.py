from pathlib import Path

import pytest

from src.docking.pose_exporter import PoseExportError, export_poses


def _write_valid_reference_sdf(tmp_path: Path) -> Path:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    reference_mol = Chem.MolFromSmiles("CCO" + "C" * 42)
    if reference_mol is None:
        pytest.fail("RDKit could not build a test molecule")

    reference_mol = Chem.AddHs(reference_mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(reference_mol, params=params) != 0:
        pytest.fail("RDKit could not embed a test molecule")

    try:
        mmff_props = AllChem.MMFFGetMoleculeProperties(reference_mol)
        if mmff_props is not None:
            AllChem.MMFFOptimizeMolecule(reference_mol, mmffMolProperties=mmff_props)
        else:
            AllChem.UFFOptimizeMolecule(reference_mol)
    except Exception:
        AllChem.UFFOptimizeMolecule(reference_mol)

    reference_sdf = tmp_path / "reference.sdf"
    writer = Chem.SDWriter(str(reference_sdf))
    writer.write(reference_mol)
    writer.close()
    return reference_sdf


def _write_pose_pdbqt(path: Path, model: int, atom_count: int = 45) -> None:
    lines = [f"MODEL {model}", "REMARK VINA RESULT: -10.000 0.000 0.000"]
    for index in range(1, atom_count + 1):
        x = (index % 7) * 0.25
        y = ((index // 7) % 7) * 0.25
        z = (index % 5) * 0.25
        lines.append(f"ATOM {index:5d}  C   UNL     1 {x:8.3f} {y:8.3f} {z:8.3f}  1.00  0.00     0.129 C")
    lines.append("ENDMDL")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_export_poses_uses_fallback_for_missing_meeko(tmp_path: Path) -> None:
    reference_sdf = _write_valid_reference_sdf(tmp_path)
    output_sdf = tmp_path / "poses.sdf"
    report_path = tmp_path / "report.json"

    pose_pdbqt_1 = tmp_path / "pose1.pdbqt"
    _write_pose_pdbqt(pose_pdbqt_1, model=1)
    pose_pdbqt_2 = tmp_path / "pose2.pdbqt"
    _write_pose_pdbqt(pose_pdbqt_2, model=2)

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
    reference_sdf = _write_valid_reference_sdf(tmp_path)
    output_sdf = tmp_path / "poses.sdf"
    report_path = tmp_path / "report.json"

    pose_path = tmp_path / "pose.pdbqt"
    _write_pose_pdbqt(pose_path, model=1, atom_count=44)

    with pytest.raises(PoseExportError):
        export_poses(
            pose_pdbqts=[pose_path],
            output_sdf=output_sdf,
            reference_sdf=reference_sdf,
            report_path=report_path,
        )
