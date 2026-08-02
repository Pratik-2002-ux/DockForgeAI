# External XAI validation

The external interface accepts CSV or Parquet structures, removes exact training overlap, flags scaffold overlap, applies the locked feature schema and model, and writes predictions plus overlap and applicability reports. It never retrains, changes the decision threshold, or overwrites locked artifacts. No independent external validation is claimed until a genuine user-supplied dataset is evaluated.

The PubChem pipeline adds an official-service discovery path for independently deposited BioAssays. Its status is reported separately because assay provenance, outcome compatibility, class balance, and chemical-domain shift determine the strength of external evidence.
