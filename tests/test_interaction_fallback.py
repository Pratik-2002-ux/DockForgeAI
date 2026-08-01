from pathlib import Path

from src.analysis.interaction_engine import analyze_complex


def test_analyze_complex_fallback_detects_catalytic_contacts(tmp_path: Path) -> None:
    complex_path = tmp_path / "complex.pdb"
    complex_path.write_text(
        """ATOM      1  N   ASP A  25      0.000  0.000  0.000  1.00  0.00           N
ATOM      2  OD1 ASP A  25      1.000  0.000  0.000  1.00  0.00           O
ATOM      3  O   LIG L  1      0.900  0.100  0.000  1.00  0.00           O
ATOM      4  C   LIG L  1      2.000  0.000  0.000  1.00  0.00           C
TER
""",
        encoding="utf-8",
    )

    result = analyze_complex(complex_path=complex_path, pose_label="crystal", pose_rank=1, engine="fallback")

    assert result["engine"] == "fallback"
    assert any(item["interaction_type"] == "hydrogen_bond" for item in result["interactions"])
    assert any(item["protein_residue_name"] == "ASP" and item["protein_residue_number"] == 25 for item in result["interactions"])
