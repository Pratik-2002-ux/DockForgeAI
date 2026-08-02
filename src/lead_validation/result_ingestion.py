def selectivity_index(efficacy,cytotoxicity,efficacy_units,cytotoxicity_units,compatible=True):
 if not compatible or efficacy_units!=cytotoxicity_units:raise ValueError("incompatible assay conditions")
 if float(efficacy)<=0:raise ValueError("efficacy must be positive")
 return float(cytotoxicity)/float(efficacy)
