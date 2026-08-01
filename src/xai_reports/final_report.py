from __future__ import annotations
import hashlib,json
from pathlib import Path
from src.xai_reports.figure_generator import generate_figures
from src.xai_reports.model_card import generate_model_cards
from src.xai_reports.html_report import generate_html
from src.xai_reports.pdf_report import generate_pdf
def generate_xai_report(root:Path):
 figures=generate_figures(root);cards=generate_model_cards(root);html=generate_html(root);pdf=generate_pdf(root);files=figures+cards+[html,pdf];manifest={"generated_files":[{"path":p.relative_to(root).as_posix(),"size_bytes":p.stat().st_size,"sha256":hashlib.sha256(p.read_bytes()).hexdigest()} for p in files],"limitations":["DUD-E decoys are not experimentally inactive","Source shortcuts can inflate performance","Explanations are not biological causality"]};path=root/"data/processed/xai_final_reports/xai_report_manifest.json";path.write_text(json.dumps(manifest,indent=2));return manifest
