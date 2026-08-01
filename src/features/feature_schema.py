"""Feature schema and leakage guards."""
FORBIDDEN={"source_database","label_source","compound_id","standardized_compound_id","scaffold_id","dataset_split","binary_label","pchembl_value","activity_median_pchembl","source_record_id","source_target_id"}
def validate_feature_names(names:list[str])->None:
    if len(names)!=len(set(names)): raise ValueError("Duplicate feature names")
    leaked=FORBIDDEN.intersection(names)
    if leaked: raise ValueError(f"Forbidden leakage features: {sorted(leaked)}")
