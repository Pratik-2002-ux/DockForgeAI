"""Common source-record schema."""

COMMON_COLUMNS = [
    "source_database", "source_record_id", "source_compound_id",
    "source_assay_id", "source_document_id", "source_publication_id",
    "source_depositor", "target_name", "target_identifier",
    "target_organism", "target_variant", "assay_type", "activity_type",
    "original_value", "original_units", "activity_relation",
    "standard_value_nm", "pactivity", "categorical_outcome",
    "original_smiles", "canonical_isomeric_smiles", "inchi", "inchikey",
    "inchikey_connectivity", "murcko_scaffold_smiles", "measurement_date",
    "publication_date", "retrieval_date", "provenance_reference",
    "source_quality_score", "assay_quality_score", "external_overlap_flag",
    "record_status", "ki_nm", "kd_nm", "ic50_nm", "ec50_nm", "pki",
    "pkd", "pic50", "pec50",
]
