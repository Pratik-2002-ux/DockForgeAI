from src.datasets.scaffold_splitter import assign_scaffold_splits, scaffold_leakage

def records(): return [{"scaffold_id":f"s{i//2}","binary_label":"active" if i%2 else "inactive_or_decoy","id":i} for i in range(20)]
def test_scaffolds_do_not_leak() -> None:
    assigned=assign_scaffold_splits(records()); assert scaffold_leakage(assigned)=={}
def test_scaffold_split_deterministic() -> None:
    assert assign_scaffold_splits(records())==assign_scaffold_splits(records())
def test_scaffold_split_approximately_tracks_fractions() -> None:
    assigned=assign_scaffold_splits(records()); counts={name:sum(row["dataset_split"]==name for row in assigned) for name in ("train","validation","test")}
    assert counts["train"] >= counts["validation"] and counts["train"] >= counts["test"]
