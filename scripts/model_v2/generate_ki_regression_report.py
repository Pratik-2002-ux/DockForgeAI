#!/usr/bin/env python3
import argparse,json
from pathlib import Path
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer
ROOT=Path(__file__).resolve().parents[2]
def main():
 p=argparse.ArgumentParser();p.add_argument("--overwrite",action="store_true");p.parse_args();data=json.loads((ROOT/"data/model_v2/regression/evaluation/best_model.json").read_text());manifest=json.loads((ROOT/"data/model_v2/regression/dataset/ki_dataset_manifest.json").read_text());out=ROOT/"reports/model_v2";out.mkdir(parents=True,exist_ok=True);html=out/"DockForgeAI_v2_Ki_Regression_Report.html";pdf=out/"DockForgeAI_v2_Ki_Regression_Report.pdf";warning="Internal scaffold evaluation on a small BindingDB-only dataset. Predictions do not establish inhibition, binding, therapeutic activity, or biological causality."
 html.write_text(f"<!doctype html><h1>DockForgeAI v2 Ki Regression</h1><p>{warning}</p><h2>Curation</h2><pre>{json.dumps(manifest,indent=2)}</pre><h2>Evaluation</h2><pre>{json.dumps(data,indent=2)}</pre>")
 s=getSampleStyleSheet();story=[Paragraph("DockForgeAI v2 Ki Regression",s["Title"]),Paragraph(warning,s["BodyText"]),Spacer(1,12),Paragraph("Curation and evaluation",s["Heading2"]),Paragraph(json.dumps({"dataset":manifest["stats"],"best_model":data["best_model"],"test":data["locked_scaffold_test_metrics"]}),s["Code"])]
 SimpleDocTemplate(str(pdf),pagesize=letter).build(story);print(html);print(pdf);return 0
if __name__=="__main__":raise SystemExit(main())
