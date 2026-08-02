"""Dataset statistics helpers."""
def readiness(summary):
    active=summary.get("final_active_count",0); inactive=summary.get("final_inactive_count",0); sources=summary.get("source_count",0); assays=summary.get("assay_count",0); scaffolds=summary.get("unique_scaffold_count",0); total=summary.get("unique_compound_count",0); overlap=summary.get("locked_external_overlap_count",0)
    classical=active>=200 and inactive>=200 and sources>=2 and assays>=10 and overlap==0
    return classical, classical and total>=500 and scaffolds>=100
