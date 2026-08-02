METHODS={"surface_binding","calorimetric_binding","thermophoretic_binding","thermal_shift","other_justified"}
def validate_method(value):
 if value not in METHODS:raise ValueError("unsupported orthogonal evidence category")
 return value
