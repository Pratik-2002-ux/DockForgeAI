from pathlib import Path

import pytest

from src.docking.rmsd import RMSDCalculationError, calculate_rmsd


def test_calculate_rmsd_merges_scores_with_rmsd(tmp_path: Path) -> None:
    reference_sdf = Path("data/processed/crystallographic_ligands/1HSG_MK1_crystallographic_corrected.sdf")
    pose_sdf = Path("data/processed/docking_poses/1HSG_MK1/1HSG_MK1_redocked_poses.sdf")
    score_table_path = Path("data/processed/docking_tables/1HSG_MK1_vina_scores.csv")
    output_table_path = tmp_path / "merged.csv"
    metadata_path = tmp_path / "metadata.json"

    metadata = calculate_rmsd(
        reference_sdf=reference_sdf,
        pose_sdf=pose_sdf,
        score_table_path=score_table_path,
        output_table_path=output_table_path,
        metadata_path=metadata_path,
        threshold=2.0,
    )

    assert metadata["pose_count"] == 10
    assert output_table_path.exists()
    assert (output_table_path.read_text(encoding="utf-8").splitlines()[1].startswith("1,"))


def test_calculate_rmsd_rejects_mismatched_pose_count(tmp_path: Path) -> None:
    reference_sdf = Path("data/processed/crystallographic_ligands/1HSG_MK1_crystallographic_corrected.sdf")
    pose_sdf = Path("data/processed/docking_poses/1HSG_MK1/1HSG_MK1_redocked_poses.sdf")
    score_table_path = tmp_path / "scores.csv"
    score_table_path.write_text("pose_rank,affinity_kcal_mol\n1,-10.0\n", encoding="utf-8")
    output_table_path = tmp_path / "merged.csv"
    metadata_path = tmp_path / "metadata.json"

    with pytest.raises(RMSDCalculationError):
        calculate_rmsd(
            reference_sdf=reference_sdf,
            pose_sdf=pose_sdf,
            score_table_path=score_table_path,
            output_table_path=output_table_path,
            metadata_path=metadata_path,
            threshold=2.0,
        )
