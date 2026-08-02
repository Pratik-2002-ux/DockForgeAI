#!/usr/bin/env python3
import json
from pathlib import Path
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer
ROOT=Path(__file__).resolve().parents[2];readiness=json.loads((ROOT/"data/model_v2/readiness/readiness_report.json").read_text());out=ROOT/"reports/model_v2";out.mkdir(parents=True,exist_ok=True);html=out/"DockForgeAI_v2_Binary_Dataset_Readiness_Report.html";pdf=out/"DockForgeAI_v2_Binary_Dataset_Readiness_Report.pdf";message="TRAINING BLOCKED" if not readiness["training_allowed"] else "DATASET READY"
html.write_text(f"<!doctype html><h1>{message}</h1><p>Binary data development only; no binary model was fitted.</p><pre>{json.dumps(readiness,indent=2)}</pre>")
s=getSampleStyleSheet();story=[Paragraph(message,s["Title"]),Paragraph("Binary data development only; no binary model was fitted.",s["BodyText"]),Spacer(1,12),Paragraph("Blockers",s["Heading2"])] + [Paragraph(x,s["BodyText"]) for x in readiness.get("blockers",[])];SimpleDocTemplate(str(pdf),pagesize=letter).build(story)
