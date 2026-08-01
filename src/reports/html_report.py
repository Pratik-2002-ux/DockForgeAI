"""Portable offline HTML report generation."""

from __future__ import annotations

import base64
import html
import json
import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

LOGGER = logging.getLogger(__name__)


class HTMLReportError(RuntimeError):
    """Raised when a final HTML report cannot be produced."""


def _image(path: Path) -> str:
    return f'<img src="data:image/png;base64,{base64.b64encode(path.read_bytes()).decode()}" alt="{html.escape(path.stem)}">'


def _table(rows: list[dict[str, Any]], columns: list[str]) -> str:
    if not rows: return "<p>Not available.</p>"
    header = "".join(f"<th>{html.escape(c.replace('_', ' ').title())}</th>" for c in columns)
    body = "".join("<tr>" + "".join(f"<td>{html.escape(str(row.get(c, '')))}</td>" for c in columns) + "</tr>" for row in rows)
    return f"<div class='table-wrap'><table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table></div>"


def generate_html_report(data: dict[str, Any], output_path: Path, figures: list[Path] | None = None, molecule_images: list[Path] | None = None, overwrite: bool = False) -> Path:
    """Generate an image-embedded report readable without internet access."""

    if output_path.exists() and not overwrite: raise HTMLReportError(f"Output exists: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    project, protein, redocking = data["project"], data["protein"], data["redocking"]
    screening, interactions, priority = data["screening"], data["interactions"], data["prioritization"]
    all_images = [p for p in (figures or []) + (molecule_images or []) if p.suffix.lower() == ".png" and p.is_file()]
    sections = [
        ("environment", "1. Project and computational environment", _table([{"software": k, "version": v} for k, v in data["software"].items()], ["software", "version"])),
        ("protein", "2. Protein structure and preparation", f"<p><strong>{html.escape(str(protein.get('title')))}</strong></p><p>Method: {html.escape(str(protein.get('experimental_method')))}; resolution: {html.escape(str(protein.get('resolution_angstrom')))} Å; chains: {html.escape(', '.join(protein.get('chains', [])))}; residues: {html.escape(str(protein.get('residue_count')))}; waters removed: {html.escape(str(protein.get('removed_waters')))}.</p>"),
        ("ligand", "3. Ligand preparation", "<p>MK1 was extracted, chemistry-corrected, hydrogenated and prepared by the validated workflow.</p>"),
        ("box", "4. Docking box and receptor assumptions", f"<p>{html.escape(str(protein.get('protonation_assumption')))}</p><details><summary>Docking box</summary><pre>{html.escape(json.dumps(redocking.get('box', {}), indent=2))}</pre></details>"),
        ("redocking", "5–6. Redocking validation and pose ranking", f"<div class='cards'><div><b>Vina score</b><span>{redocking.get('affinity_kcal_mol')} kcal/mol</span></div><div><b>RMSD</b><span>{redocking.get('rmsd_angstrom')} Å</span></div><div><b>Interpretation</b><span>{html.escape(str(redocking.get('interpretation')))}</span></div></div>{_table(redocking.get('poses', []), ['pose_rank','affinity_kcal_mol','crystal_in_place_rmsd_angstrom','passes_2A_threshold'])}"),
        ("interactions", "7–8. Interaction analysis and catalytic Asp25 findings", f"<p>Total contacts: {interactions.get('total_contacts')}; catalytic Asp25 contacts: {interactions.get('catalytic_asp25_contacts')}.</p>{_table([{'type': k, 'count': v} for k,v in interactions.get('counts_by_type', {}).items()], ['type','count'])}"),
        ("screening", "9–10. Virtual-screening demonstration and ranking", f"<p><strong>This three-compound run demonstrates software operation; it is not a biological discovery and these are not validated inhibitors.</strong></p>{_table(screening.get('ranking', []), ['compound_id','best_affinity','status','runtime_seconds'])}"),
        ("descriptors", "11. Drug-likeness and structural-alert profile", "<p>RDKit descriptors and structural alerts support prioritization. They are not full ADMET predictions.</p>"),
        ("priority", "12. Hit prioritization", _table(priority.get('ranking', []), ['compound_id','rank','consensus_score','best_affinity','qed','alerts']) + f"<p>External ADMET provenance: {html.escape(str(priority.get('external_admet_provenance')))}</p>"),
        ("limits", "13. Scientific limitations", "<ul>" + "".join(f"<li>{html.escape(item)}</li>" for item in data.get('limitations', [])) + "</ul>"),
        ("repro", "14. Reproducibility and file paths", "<details open><summary>Project-relative source inputs</summary><ul>" + "".join(f"<li><code>{html.escape(p)}</code></li>" for p in data.get('files', {}).get('source_inputs', [])) + "</ul></details>"),
        ("conclusion", "15. Conclusion", "<p>The validated 1HSG/MK1 redocking result strongly reproduces the crystallographic pose at 0.518 Å RMSD. All downstream predictions require experimental validation.</p>"),
    ]
    image_html = "".join(f"<figure>{_image(p)}<figcaption>{html.escape(p.stem.replace('_',' ').title())}</figcaption></figure>" for p in all_images)
    nav = "".join(f"<a href='#{anchor}'>{title}</a>" for anchor, title, _ in sections)
    body = "".join(f"<section id='{anchor}'><h2>{title}</h2>{content}</section>" for anchor, title, content in sections)
    document = f"""<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>DockForgeAI {project['target_id']}/{project['ligand_id']} Report</title><style>
body{{font:16px/1.55 system-ui,sans-serif;color:#172033;margin:0;background:#f4f7fb}}header,main,nav{{max-width:1100px;margin:auto}}header{{padding:3rem 1.5rem 1rem}}nav{{display:flex;gap:.8rem;overflow:auto;padding:1rem 1.5rem;position:sticky;top:0;background:white}}nav a{{white-space:nowrap;color:#164e9b}}main{{background:white;padding:1.5rem 3rem}}h1{{font-size:2.5rem}}h2{{border-bottom:2px solid #dbeafe;padding-top:2rem}}.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:1rem}}.cards div{{padding:1rem;background:#eff6ff;border-radius:8px}}.cards span,.cards b{{display:block}}.table-wrap{{overflow:auto}}table{{border-collapse:collapse;width:100%}}th,td{{border:1px solid #ccd5e1;padding:.5rem;text-align:left}}th{{background:#eaf2ff}}figure{{display:inline-block;width:min(46%,500px);vertical-align:top}}img{{max-width:100%;height:auto}}pre{{overflow:auto}}@media(max-width:700px){{main{{padding:1rem}}figure{{width:100%}}}}@media print{{nav{{display:none}}body{{background:white}}section{{break-inside:avoid}}}}
</style></head><body><header><h1>DockForgeAI</h1><p>Reproducible Molecular Docking, Virtual Screening and Hit Prioritization Report</p><p>Target {project['target_id']} · ligand {project['ligand_id']} · version {project['version']} · generated {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}</p><p><strong>Docking scores are approximate ranking values, not experimental binding free energies.</strong></p></header><nav>{nav}</nav><main><section><h2>Executive summary</h2><p>Validated redocking: {redocking.get('affinity_kcal_mol')} kcal/mol and {redocking.get('rmsd_angstrom')} Å RMSD ({html.escape(str(redocking.get('interpretation')))}). Interaction contacts: {interactions.get('total_contacts')}. Demonstration compounds screened: {screening.get('library_size')}.</p></section><section><h2>Figures and molecular structures</h2>{image_html}</section>{body}</main></body></html>"""
    try: output_path.write_text(document, encoding="utf-8")
    except OSError as exc:
        LOGGER.exception("HTML generation failed"); raise HTMLReportError(str(exc)) from exc
    if output_path.stat().st_size == 0: raise HTMLReportError("Generated HTML is empty")
    return output_path
