# PubChem external validation

DockForgeAI discovers candidate HIV-1 protease BioAssays through official NCBI Entrez and PubChem PUG services, rejects ChEMBL depositions, audits direct target measurement and outcome clarity, standardizes structures, removes exact benchmark overlap, and applies the locked model without retraining, recalibration, or threshold tuning.

The accepted BindingDB assays are biochemical literature-extracted records with compatible PubChem categorical outcomes, but their pooled evaluation remains cautionary because assay methods differ, inactive classes are sparse, and literature provenance may overlap records curated elsewhere. PubChem outcomes establish experimental assay observations, not therapeutic efficacy.

Raw API responses are locally cached and excluded from Git. See the generated discovery, dataset, evaluation, and report manifests for endpoints and checksums.
