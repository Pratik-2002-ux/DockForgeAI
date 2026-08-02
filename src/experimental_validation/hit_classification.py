STATUSES={"not tested","test failed","inconclusive","no confirmed inhibition","preliminary biochemical signal","reproducible biochemical hit","quantitative biochemical hit"}
def classify(qc_pass,replicate_pass,quantitative=False,signal=False):
 if not qc_pass:return "test failed"
 if not signal:return "no confirmed inhibition"
 if not replicate_pass:return "inconclusive"
 return "quantitative biochemical hit" if quantitative else "reproducible biochemical hit"
