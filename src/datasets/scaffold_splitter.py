"""Deterministic class-aware scaffold group splitting."""

from __future__ import annotations

import random
from collections import Counter, defaultdict
from typing import Any


def assign_scaffold_splits(records: list[dict[str, Any]], seed: int = 2026, fractions: tuple[float, float, float] = (0.70, 0.15, 0.15)) -> list[dict[str, Any]]:
    """Assign intact scaffold groups to train, validation, and test."""
    if abs(sum(fractions) - 1.0) > 1e-9: raise ValueError("Split fractions must sum to one")
    groups: dict[str, list[int]] = defaultdict(list)
    for index, row in enumerate(records): groups[row["scaffold_id"]].append(index)
    rng = random.Random(seed); items = list(groups.items()); rng.shuffle(items); items.sort(key=lambda item: len(item[1]), reverse=True)
    names = ("train", "validation", "test"); targets = {name: len(records) * fraction for name, fraction in zip(names, fractions)}
    class_targets = {name: Counter(row["binary_label"] for row in records) for name in names}
    for name, fraction in zip(names, fractions): class_targets[name] = Counter({label: count * fraction for label, count in class_targets[name].items()})
    counts = Counter(); label_counts = {name: Counter() for name in names}; assignments: dict[int, str] = {}
    for _, indices in items:
        labels = Counter(records[index]["binary_label"] for index in indices)
        def need(name: str) -> float:
            size_need = (targets[name] - counts[name]) / max(targets[name], 1)
            class_need = sum((class_targets[name][label] - label_counts[name][label]) / max(class_targets[name][label], 1) for label in labels) / max(len(labels), 1)
            return size_need + class_need
        chosen = max(names, key=lambda name: (need(name), -names.index(name)))
        for index in indices: assignments[index] = chosen
        counts[chosen] += len(indices); label_counts[chosen].update(labels)
    return [{**row, "dataset_split": assignments[index]} for index, row in enumerate(records)]


def scaffold_leakage(records: list[dict[str, Any]]) -> dict[str, list[str]]:
    """Return scaffold IDs found in multiple splits."""
    observed: dict[str, set[str]] = defaultdict(set)
    for row in records: observed[row["scaffold_id"]].add(row["dataset_split"])
    return {scaffold: sorted(splits) for scaffold, splits in observed.items() if len(splits) > 1}
