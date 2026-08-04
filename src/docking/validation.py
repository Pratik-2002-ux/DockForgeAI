"""Package existing computed crystallographic redocking validation for Phase 4."""
from __future__ import annotations
import json
from pathlib import Path
from typing import Any
from src.docking.phase4_common import atomic_json, provenance, sha256, write_checksums, write_reports


def validate_docking(summary_path: Path, analysis_path: Path, output_dir: Path, root: Path, overwrite: bool = False) -> dict[str, Any]:
    """Validate cached, computed 1HSG/MK1 redocking results without recomputation."""
    if output_dir.exists() and any(output_dir.iterdir()) and not overwrite: raise FileExistsError("Validation outputs exist; use --overwrite")
    output_dir.mkdir(parents=True, exist_ok=True); summary = json.loads(summary_path.read_text())
    if summary.get("pose_count", 0) < 1 or "min_rmsd_angstrom" not in summary: raise ValueError("Redocking summary lacks computed RMSD results")
    manifest = {"schema_version": 1, "module": "4.5_docking_validation", "redocking_summary_path": str(summary_path.relative_to(root)), "redocking_summary_checksum": sha256(summary_path), "analysis_table_path": str(analysis_path.relative_to(root)), "analysis_table_checksum": sha256(analysis_path), "pose_rmsd_angstrom": summary["rmsd_best_affinity_pose_angstrom"], "minimum_pose_rmsd_angstrom": summary["min_rmsd_angstrom"], "rmsd_threshold_angstrom": summary["threshold_angstrom"], "pose_reproduction_success": summary["overall_pose_reproduction_success"], "binding_site_coverage": "crystallographic-ligand-derived docking box", "reproducibility": "fixed Vina seed 2026; single recorded run", "warning": "A replicate-run variance estimate was not fabricated.", **provenance(root)}
    write_reports("DockForgeAI Docking Validation", {"Best-affinity pose RMSD (Å)": manifest["pose_rmsd_angstrom"], "Minimum RMSD (Å)": manifest["minimum_pose_rmsd_angstrom"], "Threshold (Å)": manifest["rmsd_threshold_angstrom"], "Pose reproduction": manifest["pose_reproduction_success"], "Binding-site coverage": manifest["binding_site_coverage"], "Reproducibility": manifest["reproducibility"]}, output_dir / "validation_report.html", output_dir / "validation_report.pdf")
    outputs = [output_dir / "validation_report.html", output_dir / "validation_report.pdf"]; manifest["output_checksums"] = {path.name: sha256(path) for path in outputs}; atomic_json(manifest, output_dir / "validation_manifest.json"); write_checksums([*outputs, output_dir / "validation_manifest.json"], output_dir / "validation_checksums.json"); return manifest
