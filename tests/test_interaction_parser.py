from pathlib import Path

from src.analysis.interaction_parser import parse_interactions


def test_parse_interactions_standardizes_records(tmp_path: Path) -> None:
    records = [
        {
            "interaction_type": "hydrogen_bond",
            "protein_chain": "A",
            "protein_residue_name": "ASP",
            "protein_residue_number": 25,
            "protein_atom": "OD1",
            "ligand_atom": "O1",
            "distance_angstrom": 2.8,
            "angle_degrees": 125.0,
            "source_engine": "fallback",
            "pose_label": "crystal",
            "pose_rank": 1,
            "notes": "predicted",
        }
    ]

    result = parse_interactions(records, output_dir=tmp_path, pose_label="crystal", pose_rank=1)

    assert result["pose_label"] == "crystal"
    assert result["row_count"] == 1
    assert (tmp_path / "crystal_interactions.csv").exists()
    assert (tmp_path / "crystal_interactions.json").exists()
