"""Strict relation-aware measurement parsing."""
from __future__ import annotations
import math, re

PATTERN = re.compile(r"^\s*(<=|>=|<|>|=)?\s*([0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)\s*$")

def parse_measurement(value: object, relation: object = None) -> tuple[str | None, float | None]:
    if value is None: return None, None
    text = str(value).strip()
    if not text or text.lower() == "nan": return None, None
    match = PATTERN.fullmatch(text)
    if not match: return None, None
    rel = str(relation or match.group(1) or "=").strip()
    number = float(match.group(2))
    if rel not in {"=", "<", "<=", ">", ">="} or not math.isfinite(number) or number <= 0:
        return None, None
    return rel, number
