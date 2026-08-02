#!/usr/bin/env python3
"""Generate honest phase-status reports without fabricating evidence."""
from pathlib import Path
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate,Paragraph
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/"reports/model_v2";OUT.mkdir(parents=True,exist_ok=True)
REPORTS={"DockForgeAI_Prospective_Screening_Report":"Phase 2 incomplete: no prospective supplier library supplied and no compounds screened.","DockForgeAI_Biochemical_Validation_Report":"Phase 3 incomplete: no laboratory study registered and no genuine biochemical results imported.","DockForgeAI_Orthogonal_Validation_Report":"Phase 4 incomplete: no genuine orthogonal evidence imported.","DockForgeAI_Lead_Validation_Report":"Phase 5 incomplete: no genuine cellular, ADME, or early-safety evidence imported.","DockForgeAI_Integrated_Discovery_Report":"Only Phase 1 computational Ki regression is complete. No candidate is an experimentally validated hit, binder, cellular hit, early lead, or drug."}
if __name__=="__main__":
 style=getSampleStyleSheet()
 for name,text in REPORTS.items():
  (OUT/f"{name}.html").write_text(f"<!doctype html><h1>{name.replace('_',' ')}</h1><p>{text}</p><p>Negative, missing, failed, and inconclusive evidence must remain visible.</p>")
  SimpleDocTemplate(str(OUT/f"{name}.pdf"),pagesize=letter).build([Paragraph(name.replace('_',' '),style["Title"]),Paragraph(text,style["BodyText"])])
