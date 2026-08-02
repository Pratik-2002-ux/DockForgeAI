"""Provenance completeness utilities."""
def missing_provenance(record): return not bool(record.get("source_record_id")) or not bool(record.get("provenance_reference"))
