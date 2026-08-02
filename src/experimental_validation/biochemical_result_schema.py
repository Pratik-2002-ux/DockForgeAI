REQUIRED=["compound_id","study_id","laboratory_id","assay_sop_id","batch_id","replicate_id","concentration","concentration_units","raw_response","normalized_response","positive_control_response","negative_control_response","vehicle_control_response","result_status","quality_control_status"]
def validate_columns(frame):
 missing=set(REQUIRED)-set(frame)
 if missing:raise ValueError(f"missing biochemical fields: {sorted(missing)}")
 return True
