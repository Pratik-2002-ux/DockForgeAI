def controls_pass(row): return row.get("quality_control_status")=="pass" and row.get("positive_control_response") is not None and row.get("negative_control_response") is not None
