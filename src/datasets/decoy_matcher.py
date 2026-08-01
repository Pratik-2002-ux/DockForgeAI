"""Deterministic physicochemical matching of DUD-E decoys to actives."""

from __future__ import annotations

import random
from typing import Any

import numpy as np

FEATURES = ("molecular_weight", "logp", "hbd", "hba", "rotatable_bonds", "tpsa", "formal_charge")


def select_matched_decoys(actives: list[dict[str, Any]], decoys: list[dict[str, Any]], count: int, seed: int = 2026) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Greedily match decoys using standardized descriptor distance."""
    if not actives or count <= 0: return [], decoys
    target = min(count, len(decoys)); combined = actives + decoys
    matrix = np.asarray([[float(row[key]) for key in FEATURES] for row in combined], dtype=float)
    scale = np.std(matrix, axis=0); scale[scale == 0] = 1.0
    center = np.mean(np.asarray([[float(row[key]) for key in FEATURES] for row in actives]), axis=0)
    rng = random.Random(seed); available = list(range(len(decoys))); rng.shuffle(available); selected_indices = []
    for index in range(target):
        active = actives[index % len(actives)]; desired = np.asarray([float(active[key]) for key in FEATURES])
        sample = available if len(available) <= 600 else rng.sample(available, 600)
        best = min(sample, key=lambda candidate: float(np.sum(np.abs((np.asarray([float(decoys[candidate][key]) for key in FEATURES]) - desired) / scale))))
        available.remove(best); selected_indices.append(best)
    selected_set = set(selected_indices)
    selected = [{**decoys[index], "matching_distance_to_active_centroid": float(np.sum(np.abs((np.asarray([float(decoys[index][key]) for key in FEATURES]) - center) / scale)))} for index in selected_indices]
    excluded = [{**row, "exclusion_reason": "not_selected_by_property_matching"} for index, row in enumerate(decoys) if index not in selected_set]
    return selected, excluded


def matching_report(actives: list[dict[str, Any]], decoys: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Compare active/decoy descriptor means and medians."""
    rows = []
    for feature in FEATURES:
        active = np.asarray([float(row[feature]) for row in actives]); decoy = np.asarray([float(row[feature]) for row in decoys])
        rows.append({"descriptor": feature, "active_mean": active.mean(), "decoy_mean": decoy.mean(), "active_median": np.median(active), "decoy_median": np.median(decoy), "standardized_mean_difference": abs(active.mean() - decoy.mean()) / (np.std(active) or 1.0)})
    return rows
