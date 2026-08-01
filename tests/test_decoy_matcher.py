from src.datasets.decoy_matcher import FEATURES, matching_report, select_matched_decoys

def row(i): return {**{key:float(i) for key in FEATURES},"id":i}
def test_selection_is_deterministic() -> None:
    active=[row(1),row(2)]; decoys=[row(i) for i in range(3,20)]
    assert [r["id"] for r in select_matched_decoys(active,decoys,5)[0]] == [r["id"] for r in select_matched_decoys(active,decoys,5)[0]]
def test_matching_report_has_all_features() -> None:
    assert len(matching_report([row(1)],[row(2)]))==len(FEATURES)
