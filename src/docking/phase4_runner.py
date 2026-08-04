"""Deterministic AutoDock Vina execution for Phase 4 review ligands."""
from __future__ import annotations

import json
import re
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any

import pandas as pd

from src.docking.phase4_common import atomic_json, provenance, sha256, write_checksums, write_reports

SCORE = re.compile(r"REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)")


class DockingError(RuntimeError):
    """Raised when a real docking engine fails or is unavailable."""


def parse_scores(path: Path) -> list[dict[str, float | int]]:
    """Parse computed Vina scores from a multi-model PDBQT."""

    matches = SCORE.findall(path.read_text(encoding="utf-8"))
    return [{"pose_rank": rank, "docking_score_kcal_mol": float(values[0]), "rmsd_lb": float(values[1]), "rmsd_ub": float(values[2])} for rank, values in enumerate(matches, 1)]


def write_best_pose(source: Path, destination: Path) -> None:
    """Extract model 1 from a Vina multi-model file."""

    lines: list[str] = []
    for line in source.read_text(encoding="utf-8").splitlines(keepends=True):
        lines.append(line)
        if line.startswith("ENDMDL"):
            break
    if not any("REMARK VINA RESULT" in line for line in lines):
        raise DockingError(f"No computed pose found in {source}")
    destination.write_text("".join(lines), encoding="utf-8")


def run_docking(ligand_dir: Path, receptor: Path, box_json: Path, output_dir: Path, root: Path, overwrite: bool = False, engine: str = "vina", seed: int = 2026) -> dict[str, Any]:
    """Dock every prepared ligand and preserve poses, logs, and provenance."""

    if engine != "vina":
        optional = shutil.which({"autodock-gpu": "autodock_gpu_128wi", "quickvina2": "qvina2"}.get(engine, engine))
        if not optional:
            raise DockingError(f"Optional docking engine unavailable: {engine}")
        raise DockingError(f"Engine adapter is not configured for {engine}; use vina")
    executable = shutil.which("vina")
    if not executable:
        raise DockingError("AutoDock Vina is unavailable; install vina=1.2.7")
    if output_dir.exists() and any(output_dir.iterdir()) and not overwrite:
        raise FileExistsError("Docking outputs exist; use --overwrite")
    if overwrite and output_dir.exists(): shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    box = json.loads(box_json.read_text())
    version = subprocess.run([executable, "--version"], capture_output=True, text=True, check=False).stdout.strip()
    rows: list[dict[str, Any]] = []
    ligand_paths = sorted(ligand_dir.glob("*.pdbqt"))
    if not ligand_paths: raise DockingError(f"No ligand PDBQT files found in {ligand_dir}")
    for index, ligand in enumerate(ligand_paths, 1):
        compound_id = ligand.stem; run_dir = output_dir / compound_id; run_dir.mkdir()
        all_poses = run_dir / "top10_poses.pdbqt"; log = run_dir / "docking_log.txt"
        command = [executable, "--receptor", str(receptor), "--ligand", str(ligand), "--center_x", str(box["center"]["x"]), "--center_y", str(box["center"]["y"]), "--center_z", str(box["center"]["z"]), "--size_x", str(box["size"]["x"]), "--size_y", str(box["size"]["y"]), "--size_z", str(box["size"]["z"]), "--exhaustiveness", "8", "--num_modes", "10", "--energy_range", "4", "--cpu", "4", "--seed", str(seed), "--out", str(all_poses)]
        print(f"Docking {index}/{len(ligand_paths)}: {compound_id}", flush=True); started = time.monotonic()
        completed = subprocess.run(command, cwd=root, capture_output=True, text=True, check=False, timeout=7200)
        log.write_text((completed.stdout or "") + "\nSTDERR\n" + (completed.stderr or ""), encoding="utf-8")
        if completed.returncode or not all_poses.is_file():
            raise DockingError(f"Vina failed for {compound_id}: {completed.stderr.strip() or completed.stdout.strip()}")
        scores = parse_scores(all_poses)
        if not scores: raise DockingError(f"Vina emitted no scores for {compound_id}")
        best = run_dir / "best_pose.pdbqt"; write_best_pose(all_poses, best)
        for score in scores:
            rows.append({"compound_id": compound_id, **score, "engine": "AutoDock Vina", "engine_version": version, "seed": seed, "runtime_seconds": round(time.monotonic() - started, 3), "best_pose_path": str(best.relative_to(root)), "top10_poses_path": str(all_poses.relative_to(root)), "log_path": str(log.relative_to(root))})
    results = pd.DataFrame(rows); results.to_csv(output_dir / "docking_results.csv", index=False)
    best_rows = results[results.pose_rank.eq(1)]
    manifest = {"schema_version": 1, "module": "4.3_docking", "engine": "AutoDock Vina", "engine_version": version, "receptor_path": str(receptor.relative_to(root)), "receptor_checksum": sha256(receptor), "ligand_directory": str(ligand_dir.relative_to(root)), "box_path": str(box_json.relative_to(root)), "box_checksum": sha256(box_json), "ligand_count": len(ligand_paths), "successful_ligand_count": best_rows.compound_id.nunique(), "pose_count": len(results), "exhaustiveness": 8, "num_modes": 10, "energy_range": 4, **provenance(root, seed)}
    write_reports("DockForgeAI Docking Report", {"Engine": manifest["engine"], "Version": version, "Ligands": len(ligand_paths), "Computed poses": len(results), "Score minimum": best_rows.docking_score_kcal_mol.min(), "Score maximum": best_rows.docking_score_kcal_mol.max()}, output_dir / "DockForgeAI_Docking_Report.html", output_dir / "DockForgeAI_Docking_Report.pdf")
    data_outputs = [output_dir / "docking_results.csv", output_dir / "DockForgeAI_Docking_Report.html", output_dir / "DockForgeAI_Docking_Report.pdf", *sorted(output_dir.glob("*/*.pdbqt")), *sorted(output_dir.glob("*/docking_log.txt"))]
    manifest["output_checksums"] = {str(path.relative_to(output_dir)): sha256(path) for path in data_outputs}
    atomic_json(manifest, output_dir / "docking_manifest.json")
    write_checksums([*data_outputs, output_dir / "docking_manifest.json"], output_dir / "docking_checksums.json")
    return manifest
