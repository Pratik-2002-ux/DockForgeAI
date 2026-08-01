def split_descriptor_importance(rows): return [row for row in rows if not row["feature"].startswith(("Morgan_","MACCS_"))]
