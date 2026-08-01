"""Generate Markdown and JSON reports for the redocking analysis."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SETTINGS_PATH = PROJECT_ROOT / "config" / "settings.yaml"


class ReportGenerationError(RuntimeError):
    """Raised when a redocking report cannot be generated."""


def load_settings(settings_path: Path = SETTINGS_PATH) -> dict[str, Any]:
    """Load project settings from YAML."""

    if not settings_path.exists():
        raise FileNotFoundError(f"Settings file not found: {settings_path}")

    with settings_path.open("r", encoding="utf-8") as handle:
        settings = yaml.safe_load(handle)

    if not isinstance(settings, dict):
        raise ValueError("Settings file must contain a YAML dictionary.")

    return settings


def generate_reports(
    analysis_csv: Path,
    summary_path: Path,
    execution_metadata: dict[str, Any],
    report_md_path: Path,
    report_json_path: Path,
    preparation_report_path: Path,
    box_report_path: Path,
    config_path: Path,
    threshold: float,
) -> dict[str, Any]:
    """Generate Markdown and JSON reports for the redocking workflow."""

    settings = load_settings()
    with analysis_csv.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    with summary_path.open("r", encoding="utf-8") as handle:
        summary = json.load(handle)
    with preparation_report_path.open("r", encoding="utf-8") as handle:
        preparation_report = json.load(handle)
    with box_report_path.open("r", encoding="utf-8") as handle:
        box_report = json.load(handle)

    lines = [
        "# MK1 redocking report",
        "",
        "## Project and software",
        f"- Project: {settings['project']['name']} {settings['project']['version']}",
        f"- Vina executable: {execution_metadata['executable']}",
        f"- Vina version: {execution_metadata['vina_version']}",
        f"- Runtime seconds: {execution_metadata['runtime_seconds']}",
        "",
        "## Inputs",
        f"- Receptor PDBQT: {preparation_report['receptor_pdbqt']}",
        f"- Ligand PDBQT: {config_path.read_text(encoding='utf-8').splitlines()[1].split('=', 1)[1].strip()}",
        f"- Docking box: {box_report_path}",
        f"- Docking config: {config_path}",
        f"- RMSD threshold: {threshold:.2f} Å",
        "",
        "## Settings",
        f"- Exhaustiveness: {settings.get('docking', {}).get('exhaustiveness', 'n/a')}",
        f"- Num modes: {settings.get('docking', {}).get('num_modes', 'n/a')}",
        f"- Energy range: {settings.get('docking', {}).get('energy_range', 'n/a')}",
        f"- CPU: {settings.get('docking', {}).get('cpu', 'n/a')}",
        f"- Seed: {settings.get('docking', {}).get('seed', 'n/a')}",
        "",
        "## Score and RMSD table",
        "",
        "| pose_rank | affinity_kcal_mol | vina_rmsd_lower_bound_angstrom | vina_rmsd_upper_bound_angstrom | crystal_in_place_rmsd_angstrom | passes_2A_threshold |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]

    for row in rows:
        lines.append(
            f"| {row['pose_rank']} | {row['affinity_kcal_mol']} | {row['vina_rmsd_lower_bound_angstrom']} | {row['vina_rmsd_upper_bound_angstrom']} | {row['crystal_in_place_rmsd_angstrom']} | {row['passes_2A_threshold']} |"
        )

    lines.extend(
        [
            "",
            "## Validation outcome",
            f"- Overall pose reproduction success: {summary['overall_pose_reproduction_success']}",
            f"- Top-ranked pose success: {summary['top_ranked_pose_success']}",
            f"- Best-energy pose rank: {summary['best_affinity_pose_rank']}",
            f"- Best-energy pose RMSD: {summary['rmsd_best_affinity_pose_angstrom']} Å",
            f"- Best-RMSD pose rank: {summary['min_rmsd_pose_rank']}",
            f"- Best-RMSD pose affinity: {summary['affinity_min_rmsd_pose_kcal_mol']} kcal/mol",
            f"- Interpretation: {summary['interpretation']}",
            "",
            "## Scientific limitations",
            "- Vina scores are approximate docking scores, not experimental free energies.",
            f"- The receptor protonation assumption was documented in {preparation_report_path}.",
            "- RMSD is computed in-place against the corrected crystallographic reference without prior minimization or alignment.",
            "",
            "## Reproducibility",
            f"- Vina command: {execution_metadata['command']}",
            f"- Analysis CSV: {analysis_csv}",
            f"- Summary JSON: {summary_path}",
        ]
    )
    report_md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    report = {
        "settings": settings,
        "execution_metadata": execution_metadata,
        "summary": summary,
        "preparation_report": preparation_report,
        "box_report": box_report,
        "rows": rows,
    }
    report_json_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report
