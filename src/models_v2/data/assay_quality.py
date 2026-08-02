"""Conservative assay-quality score."""
def assay_quality(direct=True, provenance=True): return .9 if direct and provenance else .5
