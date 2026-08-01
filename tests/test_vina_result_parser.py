from pathlib import Path

import pytest

from src.docking.result_parser import ResultParsingError, parse_vina_result


def test_parse_vina_result_creates_pose_files_and_score_table(tmp_path: Path) -> None:
    vina_output = tmp_path / "vina_output.pdbqt"
    vina_output.write_text(
        """MODEL 1
REMARK VINA RESULT: -10.500 1.200 2.300
ATOM      1  C1  LIG A   1      0.000 0.000 0.000
ENDMDL
MODEL 2
REMARK VINA RESULT: -9.500 3.400 4.500
ATOM      1  C1  LIG A   1      1.000 1.000 1.000
ENDMDL
""",
        encoding="utf-8",
    )

    metadata = parse_vina_result(
        vina_output_path=vina_output,
        output_dir=tmp_path / "poses",
        score_table_path=tmp_path / "scores.csv",
        metadata_path=tmp_path / "metadata.json",
    )

    assert metadata["poses_written"] == 2
    assert (tmp_path / "scores.csv").exists()
    assert (tmp_path / "poses" / "pose_001.pdbqt").exists()
    assert (tmp_path / "poses" / "pose_002.pdbqt").exists()


def test_parse_vina_result_rejects_missing_model_blocks(tmp_path: Path) -> None:
    vina_output = tmp_path / "vina_output.pdbqt"
    vina_output.write_text("REMARK VINA RESULT: -10.500 1.200 2.300\n", encoding="utf-8")

    with pytest.raises(ResultParsingError):
        parse_vina_result(
            vina_output_path=vina_output,
            output_dir=tmp_path / "poses",
            score_table_path=tmp_path / "scores.csv",
            metadata_path=tmp_path / "metadata.json",
        )
