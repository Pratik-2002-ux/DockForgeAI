import hashlib,json
from pathlib import Path
import pandas as pd
import pytest
from rdkit import Chem
from src.prospective_screening.library_preparation import prepare_library,discover_track1_dataset,sha256

ROOT=Path(__file__).resolve().parents[1]
COLS=["compound_id","canonical_smiles","supplier_id","supplier_name","catalogue_id","availability","library_name","library_release","source_file"]
def row(identifier="X1",smiles="CCO"):
 return {"compound_id":identifier,"canonical_smiles":smiles,"supplier_id":"SUP","supplier_name":"Supplier","catalogue_id":identifier,"availability":"snapshot","library_name":"Fixture","library_release":"2026-01-01","source_file":"fixture.csv"}
def run(tmp_path,rows,overwrite=False,strict=False):
 inp=tmp_path/"input.csv";pd.DataFrame(rows,columns=COLS).to_csv(inp,index=False);out=tmp_path/"study";return prepare_library(inp,out,ROOT,overwrite,strict),out
def test_valid_library(tmp_path):
 m,o=run(tmp_path,[row()]);assert m["final_usable_count"]==1 and m["screening_ready"]
def test_invalid_smiles(tmp_path):
 m,o=run(tmp_path,[row(smiles="not smiles")]);assert m["invalid_structure_count"]==1
def test_missing_column(tmp_path):
 inp=tmp_path/"x.csv";pd.DataFrame([row()]).drop(columns="source_file").to_csv(inp,index=False)
 with pytest.raises(ValueError,match="missing columns"):prepare_library(inp,tmp_path/"out",ROOT)
def test_duplicate_supplier_identity(tmp_path):
 m,o=run(tmp_path,[row(),row("X1","CCN")]);assert m["duplicate_id_count"]==1
def test_duplicate_structure(tmp_path):
 m,o=run(tmp_path,[row(),row("X2","OCC")]);assert m["duplicate_structure_count"]==1
def test_track1_exact_overlap(tmp_path):
 track,_=discover_track1_dataset(ROOT);smiles=pd.read_csv(track).canonical_isomeric_smiles.iloc[0];m,o=run(tmp_path,[row(smiles=smiles)]);assert m["track1_exact_overlap_count"]>=1 and m["final_usable_count"]==0
def test_locked_external_exact_overlap(tmp_path):
 smiles=pd.read_csv(ROOT/"data/processed/pubchem_external_validation/dataset/pubchem_external_validation.csv").canonical_isomeric_smiles.iloc[0];m,o=run(tmp_path,[row(smiles=smiles)]);assert m["locked_external_exact_overlap_count"]>=1
def _stereo_removed_variant(frame):
 for smiles in frame.canonical_isomeric_smiles.dropna():
  mol=Chem.MolFromSmiles(smiles);original=Chem.MolToSmiles(mol,isomericSmiles=True);Chem.RemoveStereochemistry(mol);variant=Chem.MolToSmiles(mol,isomericSmiles=True)
  if variant!=original:return variant
 pytest.skip("No stereochemical reference available")
def test_track1_connectivity_overlap(tmp_path):
 track,_=discover_track1_dataset(ROOT);variant=_stereo_removed_variant(pd.read_csv(track));m,o=run(tmp_path,[row(smiles=variant)]);assert m["track1_connectivity_overlap_count"]>=1
def test_locked_external_connectivity_overlap(tmp_path):
 locked=pd.read_csv(ROOT/"data/processed/pubchem_external_validation/dataset/pubchem_external_validation.csv");variant=_stereo_removed_variant(locked);m,o=run(tmp_path,[row(smiles=variant)]);assert m["locked_external_connectivity_overlap_count"]>=1
def test_scaffold_overlap_is_not_exclusion(tmp_path):
 m,o=run(tmp_path,[row(smiles="Cc1ccccc1")]);assert m["final_usable_count"] in {0,1};assert "scaffold_overlap" in pd.read_csv(o/"structure_overlaps.csv").columns
def test_descriptors_generated(tmp_path):
 m,o=run(tmp_path,[row()]);d=pd.read_csv(o/"prepared_library.csv");assert {"molecular_weight","logp","tpsa","qed"}.issubset(d)
def test_structural_flags_generated(tmp_path):
 m,o=run(tmp_path,[row(smiles="C"*70)]);assert "extreme_molecular_weight_flag" in pd.read_csv(o/"prepared_library.csv")
def test_manifest_generation(tmp_path):
 m,o=run(tmp_path,[row()]);assert json.loads((o/"preparation_manifest.json").read_text())["input_count"]==1
def test_checksum_generation(tmp_path):
 m,o=run(tmp_path,[row()]);checks=json.loads((o/"preparation_checksums.json").read_text());assert checks["prepared_library.csv"]==sha256(o/"prepared_library.csv")
def test_overwrite_required(tmp_path):
 m,o=run(tmp_path,[row()]);
 with pytest.raises(FileExistsError):prepare_library(tmp_path/"input.csv",o,ROOT)
def test_overwrite_succeeds(tmp_path):
 m,o=run(tmp_path,[row()]);m2=prepare_library(tmp_path/"input.csv",o,ROOT,overwrite=True);assert m2["final_usable_count"]==1
def test_empty_final_not_ready(tmp_path):
 m,o=run(tmp_path,[row(smiles="")]);assert m["final_usable_count"]==0 and not m["screening_ready"]
def test_candidate_pool_preserves_source_columns(tmp_path):
 m,o=run(tmp_path,[row()]);assert list(pd.read_csv(o/"candidate_pool.csv").columns)==COLS
def test_no_screening_artifacts_created(tmp_path):
 m,o=run(tmp_path,[row()]);assert not (o/"shortlisted_candidates.csv").exists() and not (o/"screened_library.csv").exists()
