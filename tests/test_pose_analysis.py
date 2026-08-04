"""Tests for Phase 4 pose analysis and docking validation."""
from pathlib import Path
from src.docking.analysis.interactions import _atoms
from src.docking.validation import validate_docking
ROOT=Path(__file__).resolve().parents[1]
def test_pdbqt_atom_parser()->None:
 atoms=_atoms(ROOT/"data/processed/docking_poses/1HSG_MK1/pose_001.pdbqt",ligand=True);assert len(atoms)>0
def test_validation_reports(tmp_path:Path)->None:
 result=validate_docking(ROOT/"data/processed/docking_results/1HSG_MK1_redocking_summary.json",ROOT/"data/processed/docking_tables/1HSG_MK1_redocking_analysis.csv",tmp_path,ROOT);assert result["pose_rmsd_angstrom"]==0.518;assert (tmp_path/"validation_report.pdf").exists()
