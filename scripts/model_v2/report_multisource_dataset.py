#!/usr/bin/env python3
"""Generate dataset-only HTML/PDF reports; contains no model claims."""
from __future__ import annotations
import argparse,json
from pathlib import Path
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer
from reportlab.lib.styles import getSampleStyleSheet
ROOT=Path(__file__).resolve().parents[2]
def main():
 p=argparse.ArgumentParser(); p.add_argument("--overwrite",action="store_true"); a=p.parse_args(); readiness=ROOT/"data/model_v2/readiness/readiness_report.json"
 if not readiness.exists(): raise FileNotFoundError("Run build_multisource_dataset.py first")
 data=json.loads(readiness.read_text()); out=ROOT/"reports/model_v2"; out.mkdir(parents=True,exist_ok=True); html=out/"DockForgeAI_v2_Multisource_Dataset_Report.html"; pdf=out/"DockForgeAI_v2_Multisource_Dataset_Report.pdf"
 sections=["Source overview","Target identity","Measurement types","Active/inactive definitions","Unit conversion","Structure standardization","Duplicate handling","Conflict handling","External overlap exclusion","Class balance","Source balance","Assay balance","Scaffold diversity","Split integrity","Readiness decision","Limitations","Reproducibility"]
 body="".join(f"<h2>{s}</h2><pre>{json.dumps(data,indent=2)}</pre>" if s in {"Source overview","Readiness decision"} else f"<h2>{s}</h2><p>See the manifest and exclusion tables for auditable details.</p>" for s in sections); html.write_text(f"<!doctype html><title>DockForgeAI v2 dataset report</title><h1>DockForgeAI v2 Multi-source Experimental Dataset</h1><p>This is a dataset-curation report; no model performance or therapeutic activity is claimed.</p>{body}")
 styles=getSampleStyleSheet(); story=[Paragraph("DockForgeAI v2 Multi-source Experimental Dataset",styles["Title"]),Paragraph("Dataset curation only; no model performance or therapeutic activity is claimed.",styles["BodyText"])]
 for s in sections: story += [Spacer(1,8),Paragraph(s,styles["Heading2"]),Paragraph("See machine-readable manifests and exclusion tables for details.",styles["BodyText"])]
 SimpleDocTemplate(str(pdf),pagesize=letter).build(story); print(html.relative_to(ROOT)); print(pdf.relative_to(ROOT)); return 0
if __name__=="__main__": raise SystemExit(main())
