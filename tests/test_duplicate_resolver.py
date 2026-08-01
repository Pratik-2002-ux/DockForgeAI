from src.datasets.duplicate_resolver import resolve_duplicates

def record(smiles,key,block): return {"canonical_isomeric_smiles":smiles,"inchikey":key,"inchikey_connectivity":block}
def test_duplicate_removal() -> None:
    kept,excluded=resolve_duplicates([record("CC","A","X"),record("CC","B","Y")]); assert len(kept)==1 and excluded[0]["exclusion_reason"]=="duplicate_canonical_smiles"
def test_active_overlap_prevention() -> None:
    kept,excluded=resolve_duplicates([record("CC","A","X")],{"A"}); assert not kept and excluded[0]["exclusion_reason"]=="overlaps_active_structure"
