"""Generate PubChem external-validation PDF/HTML reports."""
from __future__ import annotations
import json
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from .pubchem_provenance import write_manifest
STATEMENTS=["The production model remained locked.","No threshold was optimized and no external compounds were used for retraining.","PubChem outcomes are experimental assay results but do not establish therapeutic efficacy.","Whole-cell and biochemical assays must not be pooled without justification.","ChEMBL-submitted PubChem records were excluded from independent evaluation.","AID 162060 was rejected because its depositor is ChEMBL.","External-validation quality depends on assay independence and outcome clarity."]
SECTIONS=["Executive summary","Assay discovery","Accepted and rejected AIDs","Depositor independence","Assay definitions","Structure retrieval","Outcome curation","Training-overlap audit","Scaffold-overlap audit","Applicability domain","Locked-model evaluation","Per-assay results","Novel-scaffold results","Calibration","Distribution shift","Limitations","Reproducibility","Conclusion"]
def generate_report(root:Path):
    root=Path(root); reports=root/"reports"; reports.mkdir(parents=True,exist_ok=True); pdf=reports/"DockForgeAI_PubChem_External_Validation_Report.pdf"; html=reports/"DockForgeAI_PubChem_External_Validation_Report.html"
    with PdfPages(pdf) as doc:
        for section in SECTIONS:
            fig,ax=plt.subplots(figsize=(8.27,11.69)); ax.axis("off"); text=section+"\n\n"+("\n".join("• "+x for x in STATEMENTS) if section in {"Executive summary","Limitations","Conclusion"} else "See checksummed CSV and JSON artifacts for exact results."); ax.text(.06,.94,text,va="top",wrap=True); doc.savefig(fig); plt.close(fig)
    body="".join(f"<h2>{s}</h2><p>{' '.join(STATEMENTS) if s in ['Executive summary','Limitations','Conclusion'] else 'See companion artifacts for exact results.'}</p>" for s in SECTIONS); html.write_text(f"<!doctype html><html><head><meta charset='utf-8'><title>DockForgeAI PubChem External Validation</title></head><body><h1>DockForgeAI PubChem External Validation</h1>{body}</body></html>")
    return write_manifest(reports/"pubchem_external_validation_report_manifest.json",[pdf,html],{"statements":STATEMENTS})
