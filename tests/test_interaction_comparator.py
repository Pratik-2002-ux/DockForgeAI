from pathlib import Path

from src.analysis.interaction_comparator import compare_interactions


def test_compare_interactions_reports_shared_lost_and_gained(tmp_path: Path) -> None:
    crystal = [
        {"protein_residue_name": "ASP", "protein_residue_number": 25, "interaction_type": "hydrogen_bond"},
        {"protein_residue_name": "LYS", "protein_residue_number": 30, "interaction_type": "hydrogen_bond"},
    ]
    docked = [
        {"protein_residue_name": "ASP", "protein_residue_number": 25, "interaction_type": "hydrogen_bond"},
        {"protein_residue_name": "GLU", "protein_residue_number": 40, "interaction_type": "hydrophobic"},
    ]

    result = compare_interactions(crystal=crystal, docked=docked, output_path=tmp_path / "comparison.json")

    assert result["shared_residue_contacts"] == ["A:25", "A:30"]
    assert result["lost_residue_contacts"] == ["A:30"]
    assert result["gained_residue_contacts"] == ["A:40"]
    assert result["jaccard_similarity"] > 0.0
