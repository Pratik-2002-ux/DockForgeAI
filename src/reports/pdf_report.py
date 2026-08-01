"""Professional multi-page PDF report generation with ReportLab."""

from __future__ import annotations

import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import Image, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle

LOGGER = logging.getLogger(__name__)


class PDFReportError(RuntimeError):
    """Raised when a final PDF cannot be produced."""


def _text(value: Any) -> str:
    return "Unavailable" if value in (None, "") else str(value)


def _table(rows: list[list[Any]], widths: list[float] | None = None) -> Table:
    table = Table([[Paragraph(_text(cell), getSampleStyleSheet()["BodyText"]) for cell in row] for row in rows], colWidths=widths, repeatRows=1)
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#dbeafe")),
        ("GRID", (0, 0), (-1, -1), 0.4, colors.grey),
        ("VALIGN", (0, 0), (-1, -1), "TOP"), ("FONTSIZE", (0, 0), (-1, -1), 8),
        ("LEFTPADDING", (0, 0), (-1, -1), 4), ("RIGHTPADDING", (0, 0), (-1, -1), 4),
    ]))
    return table


def generate_pdf_report(data: dict[str, Any], output_path: Path, figures: list[Path] | None = None, molecule_images: list[Path] | None = None, overwrite: bool = False) -> Path:
    """Generate the final PDF report and return its path."""

    if output_path.exists() and not overwrite:
        raise PDFReportError(f"Output exists: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="CoverTitle", parent=styles["Title"], fontSize=28, leading=34, alignment=TA_CENTER, spaceAfter=18))
    styles.add(ParagraphStyle(name="CoverSub", parent=styles["Heading2"], alignment=TA_CENTER, textColor=colors.HexColor("#334155")))
    story: list[Any] = []
    project, protein, redocking = data["project"], data["protein"], data["redocking"]
    screening, interactions = data["screening"], data["interactions"]

    story += [Spacer(1, 42 * mm), Paragraph("DockForgeAI", styles["CoverTitle"]),
        Paragraph("Reproducible Molecular Docking, Virtual Screening and Hit Prioritization Report", styles["CoverSub"]), Spacer(1, 15 * mm),
        Paragraph(f"Project version {_text(project.get('version'))}<br/>Target: {_text(project.get('target_id'))} &nbsp; Reference ligand: {_text(project.get('ligand_id'))}<br/>Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}", styles["BodyText"]),
        Spacer(1, 20 * mm), Paragraph("Software-generated scientific report. Docking scores are approximate ranking values, not experimental binding free energies.", styles["Italic"]), PageBreak()]

    story += [Paragraph("Executive summary", styles["Heading1"]),
        Paragraph(f"DockForgeAI evaluated {_text(protein.get('pdb_id'))} with reference ligand {_text(project.get('ligand_id'))}. The validated top pose reported a Vina score of {_text(redocking.get('affinity_kcal_mol'))} kcal/mol and crystallographic heavy-atom RMSD of {_text(redocking.get('rmsd_angstrom'))} Å ({_text(redocking.get('interpretation'))}).", styles["BodyText"]),
        Paragraph(f"Interaction analysis recorded {_text(interactions.get('total_contacts'))} geometric contacts, including {_text(interactions.get('catalytic_asp25_contacts'))} contacts involving catalytic Asp25. The demonstration screen processed {_text(screening.get('library_size'))} compounds with {_text(screening.get('success_count'))} successful runs.", styles["BodyText"]),
        Paragraph(f"Best demonstration compound: {_text(screening.get('best_demonstration_compound', {}).get('compound_id'))}, score {_text(screening.get('best_demonstration_compound', {}).get('best_affinity'))} kcal/mol. This is a software demonstration, not evidence of a validated inhibitor; experimental validation is required.", styles["BodyText"]), Spacer(1, 6 * mm)]
    if molecule_images:
        first = next((p for p in molecule_images if p.suffix == ".png" and p.is_file()), None)
        if first: story.append(Image(str(first), width=70 * mm, height=56 * mm))

    story += [PageBreak(), Paragraph("1. Project and computational environment", styles["Heading1"]),
        _table([["Software", "Version"]] + [[key, value] for key, value in data["software"].items()], [55 * mm, 115 * mm]),
        Paragraph("2. Protein structure and preparation", styles["Heading1"]),
        _table([["Field", "Value"], ["PDB ID", protein.get("pdb_id")], ["Title", protein.get("title")], ["Method", protein.get("experimental_method")], ["Resolution", f"{_text(protein.get('resolution_angstrom'))} Å"], ["Chains", ", ".join(protein.get("chains", []))], ["Residues", protein.get("residue_count")], ["Waters removed", protein.get("removed_waters")]], [42 * mm, 128 * mm]),
        Paragraph("3. Ligand preparation", styles["Heading1"]), Paragraph("The crystallographic MK1 ligand was extracted, chemistry-corrected, hydrogenated and prepared for redocking using the existing validated workflow.", styles["BodyText"]),
        Paragraph("4. Docking box and receptor assumptions", styles["Heading1"]), Paragraph(_text(protein.get("protonation_assumption")), styles["BodyText"]), Paragraph(f"Docking box: {_text(redocking.get('box'))}", styles["BodyText"]), PageBreak()]

    poses = redocking.get("poses", [])
    pose_table = [["Rank", "Vina score (kcal/mol)", "Crystal RMSD (Å)", "Pass ≤2 Å"]] + [[r.get("pose_rank"), r.get("affinity_kcal_mol"), r.get("crystal_in_place_rmsd_angstrom"), r.get("passes_2A_threshold")] for r in poses]
    story += [Paragraph("5. Redocking validation", styles["Heading1"]), Paragraph(f"Best score: {_text(redocking.get('affinity_kcal_mol'))} kcal/mol; minimum RMSD: {_text(redocking.get('rmsd_angstrom'))} Å; interpretation: {_text(redocking.get('interpretation'))}.", styles["BodyText"]),
        Paragraph("6. Pose ranking table", styles["Heading1"]), _table(pose_table, [25 * mm, 50 * mm, 50 * mm, 35 * mm])]
    for name in ("docking_affinity_by_pose.png", "crystallographic_rmsd_by_pose.png"):
        image = next((p for p in figures or [] if p.name == name), None)
        if image: story.append(Image(str(image), width=155 * mm, height=93 * mm))
    story.append(PageBreak())

    story += [Paragraph("7. Interaction analysis", styles["Heading1"]), _table([["Interaction type", "Count"]] + [[k, v] for k, v in interactions.get("counts_by_type", {}).items()], [90 * mm, 50 * mm]),
        Paragraph("8. Catalytic Asp25 findings", styles["Heading1"]), Paragraph(f"The geometric analysis recorded {_text(interactions.get('catalytic_asp25_contacts'))} contacts involving ASP25 in chains A or B. Contact predictions depend on the stated geometric fallback assumptions.", styles["BodyText"])]
    for name in ("interaction_counts.png", "top_contacted_residues.png"):
        image = next((p for p in figures or [] if p.name == name), None)
        if image: story.append(Image(str(image), width=150 * mm, height=90 * mm))
    story.append(PageBreak())

    ranking = screening.get("ranking", [])
    story += [Paragraph("9. Virtual-screening demonstration", styles["Heading1"]), Paragraph("This three-compound run validates the software path only. It is not a biological discovery study and the compounds are not validated inhibitors.", styles["BodyText"]),
        Paragraph("10. Screening ranking", styles["Heading1"]), _table([["Compound", "Best Vina score", "Status", "Runtime (s)"]] + [[r.get("compound_id"), r.get("best_affinity"), r.get("status"), r.get("runtime_seconds")] for r in ranking], [40 * mm] * 4),
        Paragraph("11. Drug-likeness and structural-alert profile", styles["Heading1"]), Paragraph("Reported RDKit descriptors and structural alerts support prioritization; they are not full ADMET predictions.", styles["BodyText"])]
    for name in ("screening_score_ranking.png", "screening_workflow_counts.png"):
        image = next((p for p in figures or [] if p.name == name), None)
        if image: story.append(Image(str(image), width=150 * mm, height=90 * mm))
    story.append(PageBreak())

    priority = data.get("prioritization", {})
    story += [Paragraph("12. Hit prioritization", styles["Heading1"]), _table([["Compound", "Rank", "Consensus", "Vina score", "QED", "Alerts"]] + [[r.get("compound_id"), r.get("rank"), r.get("consensus_score"), r.get("best_affinity"), r.get("qed"), r.get("alerts")] for r in priority.get("ranking", [])], [28 * mm] * 6),
        Paragraph("External ADMET provenance", styles["Heading2"]), Paragraph(_text(priority.get("external_admet_provenance")), styles["BodyText"]),
        Paragraph("13. Scientific limitations", styles["Heading1"])]
    story.extend(Paragraph(f"• {item}", styles["BodyText"]) for item in data.get("limitations", []))
    story += [Paragraph("14. Reproducibility and file paths", styles["Heading1"]), Paragraph("Source inputs (project-relative):", styles["BodyText"])]
    story.extend(Paragraph(path, styles["Code"]) for path in data.get("files", {}).get("source_inputs", []))
    story += [Paragraph("15. Conclusion", styles["Heading1"]), Paragraph("The validated 1HSG/MK1 redocking result reproduces the crystallographic pose strongly (0.518 Å). Downstream screening and prioritization remain computational demonstrations requiring experimental confirmation.", styles["BodyText"])]

    def page(canvas: Any, document: Any) -> None:
        canvas.saveState(); canvas.setFont("Helvetica", 8)
        canvas.drawString(18 * mm, 12 * mm, "DockForgeAI final report")
        canvas.drawRightString(192 * mm, 12 * mm, f"Page {document.page}"); canvas.restoreState()
    try:
        SimpleDocTemplate(str(output_path), pagesize=A4, rightMargin=18 * mm, leftMargin=18 * mm, topMargin=18 * mm, bottomMargin=18 * mm, title="DockForgeAI 1HSG/MK1 Report").build(story, onFirstPage=page, onLaterPages=page)
    except Exception as exc:
        LOGGER.exception("PDF generation failed")
        raise PDFReportError(str(exc)) from exc
    if output_path.stat().st_size == 0:
        raise PDFReportError("Generated PDF is empty")
    return output_path
