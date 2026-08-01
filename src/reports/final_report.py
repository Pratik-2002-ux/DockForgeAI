"""Orchestrate final DockForgeAI PDF/HTML report generation."""

from __future__ import annotations

import json
import logging
from datetime import datetime, timezone
from pathlib import Path
from time import monotonic
from typing import Any

import yaml

from src.reports.html_report import generate_html_report
from src.reports.pdf_report import generate_pdf_report
from src.reports.report_data import PROJECT_ROOT, collect_report_data
from src.reports.report_figures import generate_report_figures
from src.visualization.molecule_renderer import MoleculeRenderError, molecule_from_sdf, molecule_from_smiles, save_molecule_images

LOGGER = logging.getLogger(__name__)


class FinalReportError(RuntimeError):
    """Raised when final-report orchestration fails."""


def _configured_paths(root: Path) -> dict[str, Path]:
    settings = yaml.safe_load((root / "config" / "settings.yaml").read_text(encoding="utf-8")) or {}
    directories = settings.get("directories", {})
    return {key: root / directories[key] for key in ("report_figures", "molecule_images", "pdf_reports", "html_reports", "final_reports")}


def _render_molecules(data: dict[str, Any], root: Path, output_dir: Path, warnings: list[str]) -> list[Path]:
    outputs: list[Path] = []
    candidates: list[tuple[str, Any, str | None]] = []
    try:
        mk1 = molecule_from_sdf(root / "data/processed/crystallographic_ligands/1HSG_MK1_crystallographic_corrected.sdf")
        candidates.append(("MK1_reference", mk1, "crystallographic reference"))
    except MoleculeRenderError as exc: warnings.append(str(exc))
    for row in data["screening"].get("ranking", []):
        compound = str(row.get("compound_id", "unknown"))
        try:
            molecule = molecule_from_sdf(root / f"data/processed/screening_library/sdf/{compound}.sdf")
            candidates.append((f"compound_{compound}", molecule, f"Vina {row.get('best_affinity')} kcal/mol"))
        except MoleculeRenderError as exc: warnings.append(str(exc))
    seen: set[str] = set()
    for index, row in enumerate(data["prioritization"].get("ranking", []), start=1):
        name = f"prioritized_{row.get('compound_id', index)}_{index}"
        if not row.get("smiles") or name in seen: continue
        try: candidates.append((name, molecule_from_smiles(str(row["smiles"])), f"rank {row.get('rank', index)}"))
        except MoleculeRenderError as exc: warnings.append(str(exc))
    for name, molecule, label in candidates:
        seen.add(name)
        try: outputs.extend(save_molecule_images(molecule, output_dir / name, name, label))
        except (MoleculeRenderError, OSError) as exc:
            LOGGER.warning("Molecule rendering failed for %s: %s", name, exc); warnings.append(str(exc))
    return outputs


def generate_final_reports(
    project_root: Path = PROJECT_ROOT,
    target_id: str = "1HSG",
    ligand_id: str = "MK1",
    output_prefix: str | None = None,
    generate_pdf: bool = True,
    generate_html: bool = True,
    overwrite: bool = False,
    strict: bool = False,
    tests_status: str | None = None,
) -> dict[str, Any]:
    """Generate assets, reports and a machine-readable manifest."""

    started = monotonic(); root = project_root.resolve(); paths = _configured_paths(root)
    data = collect_report_data(root, target_id=target_id, ligand_id=ligand_id)
    warnings = list(data.get("warnings", []))
    if strict and warnings: raise FinalReportError("Strict mode rejected missing or invalid optional inputs: " + "; ".join(warnings))
    prefix = output_prefix or f"DockForgeAI_{target_id}_{ligand_id}_report"
    figures = generate_report_figures(data, paths["report_figures"])
    molecules = _render_molecules(data, root, paths["molecule_images"], warnings)
    generated: list[Path] = figures + molecules
    if generate_pdf: generated.append(generate_pdf_report(data, paths["pdf_reports"] / f"{prefix}.pdf", figures, molecules, overwrite))
    if generate_html: generated.append(generate_html_report(data, paths["html_reports"] / f"{prefix}.html", figures, molecules, overwrite))
    paths["final_reports"].mkdir(parents=True, exist_ok=True)
    manifest_path = paths["final_reports"] / "report_manifest.json"
    manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "runtime_seconds": round(monotonic() - started, 3),
        "generated_files": [{"path": p.relative_to(root).as_posix(), "size_bytes": p.stat().st_size, "timestamp": datetime.fromtimestamp(p.stat().st_mtime, timezone.utc).isoformat()} for p in generated if p.is_file()],
        "software_versions": data["software"], "source_inputs": data["files"]["source_inputs"],
        "warnings": warnings, "missing_optional_sections": data.get("warnings", []),
        "git_commit_hash": data["software"].get("git_commit"),
    }
    if tests_status is not None: manifest["tests_status"] = tests_status
    if manifest_path.exists() and not overwrite: raise FinalReportError(f"Output exists: {manifest_path}")
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    manifest["manifest_path"] = manifest_path.relative_to(root).as_posix()
    return manifest
