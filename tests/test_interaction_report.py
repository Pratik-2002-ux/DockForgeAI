from pathlib import Path

from src.analysis.interaction_report import generate_interaction_report


def test_generate_interaction_report_writes_outputs(tmp_path: Path) -> None:
    rows = [
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

    report = generate_interaction_report(rows=rows, output_dir=tmp_path)

    assert (tmp_path / "report.md").exists()
    assert (tmp_path / "report.json").exists()
    assert (tmp_path / "interactions.csv").exists()
