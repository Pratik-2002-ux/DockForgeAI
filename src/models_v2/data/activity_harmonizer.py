"""Endpoint-preserving unit conversion and binary policy."""
from __future__ import annotations
import math

UNIT_FACTORS = {"pm": .001, "nm": 1., "um": 1000., "µm": 1000., "μm": 1000., "mm": 1_000_000., "m": 1_000_000_000.}
SUPPORTED_ENDPOINTS = {"KI", "KD", "IC50", "EC50"}

def convert_to_nm(value: float, units: str) -> float:
    key = str(units).strip().lower()
    if key not in UNIT_FACTORS: raise ValueError("unsupported_units")
    result = float(value) * UNIT_FACTORS[key]
    if not math.isfinite(result) or result <= 0: raise ValueError("invalid_measurement")
    return result

def classify_measurement(value_nm: float | None, relation: str | None = "=", categorical: str | None = None, endpoint: str = "KI") -> str:
    outcome = str(categorical or "").strip().lower()
    if outcome in {"active", "inactive"}: return outcome
    if endpoint.upper() not in SUPPORTED_ENDPOINTS or value_nm is None: return "ambiguous"
    if relation in {"=", "<", "<="} and value_nm <= 1000: return "active"
    if relation in {"=", ">", ">="} and value_nm >= 10000: return "inactive"
    return "ambiguous"

def endpoint_fields(endpoint: str, value_nm: float) -> dict[str, float | None]:
    result = {k: None for k in ("ki_nm", "kd_nm", "ic50_nm", "ec50_nm", "pki", "pkd", "pic50", "pec50")}
    key = endpoint.lower(); result[f"{key}_nm"] = value_nm; result[f"p{key}"] = 9 - math.log10(value_nm)
    return result
