"""Curate and aggregate ChEMBL HIV-1 protease activity records."""

from __future__ import annotations

import math
import statistics
from collections import defaultdict
from typing import Any

from src.datasets.structure_standardizer import StructureStandardizationError, standardize_structure

UNIT_TO_NM = {"nm": 1.0, "um": 1000.0, "µm": 1000.0, "μm": 1000.0, "mm": 1_000_000.0, "m": 1_000_000_000.0, "pm": 0.001}


def _validity_reason(comment: Any) -> str | None:
    text = str(comment or "").strip().lower()
    if not text: return None
    if "outside typical range" in text: return "outside_typical_range"
    if "ambiguous" in text: return "ambiguous_measurement_annotation"
    if "conflict" in text: return "conflicting_assay_annotation"
    if "invalid" in text or "transcription error" in text: return "invalid_measurement"
    return "chembl_data_validity_flag"


def _numeric(value: Any) -> float | None:
    try:
        number = float(value)
        return number if math.isfinite(number) and number > 0 else None
    except (TypeError, ValueError): return None


def _potency(record: dict[str, Any]) -> tuple[float | None, float | None]:
    pchembl = _numeric(record.get("pchembl_value"))
    value = _numeric(record.get("standard_value"))
    unit = str(record.get("standard_units") or "").strip().lower()
    nm = value * UNIT_TO_NM[unit] if value is not None and unit in UNIT_TO_NM else None
    if pchembl is None and nm is not None: pchembl = 9.0 - math.log10(nm)
    return pchembl, nm


def curate_activities(records: list[dict[str, Any]], config: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    """Apply target, assay, potency, ambiguity, and conflict rules."""
    target = config.get("target_chembl_id", "CHEMBL243"); types = set(config.get("standard_activity_types", ["IC50", "Ki", "Kd", "EC50"]))
    active_threshold = float(config.get("active_pchembl_threshold", 6.0)); ambiguous_lower = float(config.get("ambiguous_lower_threshold", 5.0))
    minimum_confidence = int(config.get("minimum_target_confidence_score", 8))
    valid, invalid, ambiguous = [], [], []
    for index, raw in enumerate(records):
        record = dict(raw); record["source_record_id"] = str(record.get("activity_id") or f"chembl-row-{index + 1}")
        reason = None
        if record.get("target_chembl_id") != target: reason = "wrong_target"
        elif _validity_reason(record.get("data_validity_comment")): reason = _validity_reason(record.get("data_validity_comment"))
        elif record.get("standard_type") not in types: reason = "unsupported_activity_type"
        elif str(record.get("standard_relation") or "") != "=": reason = "non_exact_relation"
        elif record.get("standard_flag") not in (None, "", 1, "1"): reason = "not_standardized_activity"
        elif record.get("confidence_score") not in (None, "") and int(record["confidence_score"]) < minimum_confidence: reason = "low_target_confidence"
        pchembl, nm = _potency(record)
        smiles = record.get("canonical_smiles") or record.get("molecule_structures", {}).get("canonical_smiles")
        if reason is None and (pchembl is None or smiles is None): reason = "missing_numeric_activity_or_structure"
        if reason:
            record["exclusion_reason"] = reason; invalid.append(record); continue
        comment = str(record.get("activity_comment") or "").lower()
        if pchembl >= active_threshold: label = "active"
        elif pchembl > ambiguous_lower: label = "ambiguous"
        else: label = "inactive"
        if (label == "active" and "inactive" in comment) or (label == "inactive" and "active" in comment and "inactive" not in comment):
            record["exclusion_reason"] = "conflicting_assay_annotation"; invalid.append(record); continue
        try: structure = standardize_structure(str(smiles))
        except StructureStandardizationError as exc:
            record["exclusion_reason"] = str(exc); invalid.append(record); continue
        curated = {**record, **structure, "original_smiles": smiles, "pchembl_numeric": pchembl, "standard_activity_value_nm": nm, "record_label": label}
        if label == "ambiguous": ambiguous.append(curated)
        else: valid.append(curated)

    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in valid: grouped[record["inchikey"]].append(record)
    active, inactive, conflicting, duplicates = [], [], [], []
    for inchikey, group in grouped.items():
        labels = {row["record_label"] for row in group}; pvalues = [row["pchembl_numeric"] for row in group]
        spread = max(pvalues) - min(pvalues); representative = sorted(group, key=lambda row: row["source_record_id"])[0]
        aggregate = {**representative, "activity_record_count": len(group), "activity_min_pchembl": min(pvalues), "activity_max_pchembl": max(pvalues), "activity_median_pchembl": statistics.median(pvalues), "activity_spread": spread, "activity_conflict_flag": len(labels) > 1 or spread > 2.0}
        if len(group) > 1: duplicates.extend({"inchikey": inchikey, "source_record_id": row["source_record_id"], "record_count": len(group)} for row in group[1:])
        if aggregate["activity_conflict_flag"]: aggregate["exclusion_reason"] = "active_inactive_conflict" if len(labels) > 1 else "pchembl_range_exceeds_2"; conflicting.append(aggregate)
        elif labels == {"active"}: active.append(aggregate)
        else: inactive.append(aggregate)
    return {"active": active, "inactive": inactive, "valid_records": valid, "invalid": invalid, "ambiguous": ambiguous, "conflicting": conflicting, "duplicates": duplicates}
