from __future__ import annotations
from pathlib import Path
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate,Paragraph,PageBreak
def generate_pdf(root:Path)->Path:
 out=root/"data/processed/xai_final_reports/DockForgeAI_XAI_Benchmark_Report.pdf";out.parent.mkdir(parents=True,exist_ok=True);styles=getSampleStyleSheet();story=[]
 for i,title in enumerate(["Executive summary","Dataset and scaffold splitting","Features and models","Calibration and applicability domain","Global and local explanations","Benchmarking and limitations","Reproducibility and conclusion"]):
  story.extend([Paragraph(title,styles["Heading1"]),Paragraph("DUD-E decoys are computational decoys, not experimentally verified inactive molecules. Model explanations describe model behaviour, not biological causality. Good classification does not prove binding or therapeutic activity.",styles["BodyText"])]);story.append(PageBreak())
 SimpleDocTemplate(str(out),pagesize=A4,title="DockForgeAI XAI Benchmark").build(story);return out
