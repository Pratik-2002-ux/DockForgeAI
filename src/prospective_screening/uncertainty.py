"""Use only uncertainty information persisted with the frozen model."""
from __future__ import annotations

import numpy as np


def uncertainty_columns(predictions: np.ndarray, metadata: dict) -> tuple[dict, dict]:
    uncertainty = metadata.get("best_model", {}).get("uncertainty", {})
    half_width = uncertainty.get("conformal_half_width_95")
    if half_width is None:
        values = np.full(len(predictions), np.nan)
        return {
            "prediction_uncertainty": values,
            "prediction_interval_lower_pki": values,
            "prediction_interval_upper_pki": values,
            "prediction_interval_lower_ki_nm": values,
            "prediction_interval_upper_ki_nm": values,
            "uncertainty_status": np.repeat("unavailable", len(predictions)),
        }, {"method": "unavailable", "available": False}
    lower = predictions - float(half_width)
    upper = predictions + float(half_width)
    return {
        "prediction_uncertainty": np.full(len(predictions), np.nan),
        "prediction_interval_lower_pki": lower,
        "prediction_interval_upper_pki": upper,
        "prediction_interval_lower_ki_nm": 10 ** (9 - upper),
        "prediction_interval_upper_ki_nm": 10 ** (9 - lower),
        "uncertainty_status": np.repeat("conformal_interval_available;ensemble_std_unavailable", len(predictions)),
    }, {"method": "conformal_interval_95", "available": True, "half_width_pki": float(half_width), "per_compound_ensemble_std": False}
