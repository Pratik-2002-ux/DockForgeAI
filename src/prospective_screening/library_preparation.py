"""Deterministic prospective-library preparation without model inference."""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

from src.datasets.structure_standardizer import StructureStandardizationError, standardize_structure
from .library_validator import validate_library


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", suffix=".csv", dir=path.parent, delete=False, encoding="utf-8", newline="") as handle:
        temp = Path(handle.name); frame.to_csv(handle, index=False)
    os.replace(temp, path)


def atomic_json(value: object, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", suffix=".json", dir=path.parent, delete=False, encoding="utf-8") as handle:
        temp = Path(handle.name); json.dump(value, handle, indent=2); handle.write("\n")
    os.replace(temp, path)


def atomic_text(value: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", suffix=path.suffix, dir=path.parent, delete=False, encoding="utf-8") as handle:
        temp = Path(handle.name)
        handle.write(value)
    os.replace(temp, path)


def discover_track1_dataset(root: Path) -> tuple[Path, dict]:
    manifests = [Path(line) for line in subprocess.run(["git", "ls-files", "data/model_v2/regression/**/ki_dataset_manifest.json"], cwd=root, text=True, capture_output=True, check=True).stdout.splitlines()]
    if len(manifests) != 1: raise FileNotFoundError(f"Expected one committed Track 1 manifest, found {len(manifests)}")
    manifest_path = root / manifests[0]; manifest = json.loads(manifest_path.read_text())
    if manifest.get("endpoint") != "Ki" or manifest.get("ordinary_regression_relations") != ["="]: raise ValueError("Track 1 manifest is not authoritative exact-Ki data")
    checks_path = manifest_path.with_name("ki_dataset_checksums.json")
    if not checks_path.is_file(): raise FileNotFoundError(f"Track 1 checksums unavailable: {checks_path}")
    checks = json.loads(checks_path.read_text()); candidates = [root / key for key in checks if key.endswith("/ki_regression_dataset.csv")]
    if len(candidates) != 1 or not candidates[0].is_file(): raise FileNotFoundError("Track 1 dataset cannot be resolved from its checksum manifest")
    dataset = candidates[0]; expected = checks[str(dataset.relative_to(root))]
    if sha256(dataset) != expected: raise ValueError("Track 1 dataset checksum mismatch")
    return dataset, manifest


def _reference_maps(frame: pd.DataFrame) -> dict[str, dict[str, str]]:
    ids = frame.get("compound_id", pd.Series(frame.index.astype(str), index=frame.index)).astype(str)
    return {column: dict(zip(frame[column].dropna().astype(str), ids.loc[frame[column].notna()])) if column in frame else {} for column in ("canonical_isomeric_smiles", "inchikey", "inchikey_connectivity", "murcko_scaffold_smiles")}


def _unresolved_mixture(smiles: str) -> bool:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: return False
    organic = [fragment for fragment in Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True) if any(atom.GetAtomicNum() == 6 for atom in fragment.GetAtoms())]
    return len(organic) > 1


def prepare_library(input_path: Path, output_dir: Path, root: Path, overwrite: bool = False, strict: bool = False) -> dict:
    # RDKit emits one line per normalization warning for some supplier records;
    # keep the explicit progress output readable while retaining failures in the
    # exclusion report.
    RDLogger.DisableLog("rdApp.*")
    outputs = [output_dir / name for name in ("candidate_pool.csv", "prepared_library.csv", "excluded_compounds.csv", "structure_overlaps.csv", "descriptor_summary.csv", "preparation_manifest.json", "preparation_schema.json", "preparation_checksums.json", "preparation_summary.md")]
    existing = [path for path in outputs if path.exists()]
    if existing and not overwrite: raise FileExistsError(f"Outputs already exist; use --overwrite: {existing[0]}")
    source = pd.read_csv(input_path, dtype=str, keep_default_na=False); validate_library(source)
    track_path, _ = discover_track1_dataset(root); track = pd.read_csv(track_path); locked_path = root / "data/processed/pubchem_external_validation/dataset/pubchem_external_validation.csv"
    if not locked_path.is_file(): raise FileNotFoundError(f"Locked external dataset unavailable: {locked_path}")
    references = {"track1_ki_regression": _reference_maps(track), "locked_external": _reference_maps(pd.read_csv(locked_path))}
    output_dir.mkdir(parents=True, exist_ok=True); atomic_csv(source, output_dir / "candidate_pool.csv")
    duplicate_id = source.duplicated("compound_id", keep="first") | source.duplicated("catalogue_id", keep="first")
    rows=[]; exclusions=[]; overlaps=[]; seen_structures={}
    for position, (_, record) in enumerate(source.iterrows(), 1):
        if position == 1 or position % 1000 == 0 or position == len(source): print(f"Preparing {position:,}/{len(source):,}", flush=True)
        base={key:record[key] for key in source.columns}; original=record["canonical_smiles"]
        def exclude(reason, detail): exclusions.append({"compound_id":record["compound_id"],"original_smiles":original,"exclusion_reason":reason,"exclusion_detail":detail,"supplier_id":record["supplier_id"],"catalogue_id":record["catalogue_id"]})
        if duplicate_id.iloc[position-1]: exclude("duplicate_supplier_id", "Repeated compound_id or catalogue_id"); continue
        if not original.strip(): exclude("invalid_or_empty_smiles", "SMILES is empty"); continue
        if _unresolved_mixture(original): exclude("unresolved_mixture", "More than one organic fragment"); continue
        try: standardized=standardize_structure(original)
        except StructureStandardizationError as exc: exclude("invalid_or_empty_smiles", str(exc)); continue
        canonical=standardized["canonical_isomeric_smiles"]
        if canonical in seen_structures: exclude("duplicate_standardized_structure", f"Duplicates {seen_structures[canonical]}"); continue
        seen_structures[canonical]=record["compound_id"]
        exclusion_overlap=None
        for reference_name,maps in references.items():
            exact_id=maps["canonical_isomeric_smiles"].get(canonical) or maps["inchikey"].get(standardized["inchikey"]); connectivity_id=maps["inchikey_connectivity"].get(standardized["inchikey_connectivity"]); scaffold_id=maps["murcko_scaffold_smiles"].get(standardized["murcko_scaffold_smiles"])
            if exact_id: overlaps.append({"compound_id":record["compound_id"],"canonical_isomeric_smiles":canonical,"overlap_type":"exact","reference_dataset":reference_name,"reference_compound_id":exact_id,"exact_overlap":True,"connectivity_overlap":True,"scaffold_overlap":bool(scaffold_id)}); exclusion_overlap=exclusion_overlap or ("exact_reference_overlap",reference_name)
            elif connectivity_id: overlaps.append({"compound_id":record["compound_id"],"canonical_isomeric_smiles":canonical,"overlap_type":"connectivity","reference_dataset":reference_name,"reference_compound_id":connectivity_id,"exact_overlap":False,"connectivity_overlap":True,"scaffold_overlap":bool(scaffold_id)}); exclusion_overlap=exclusion_overlap or ("connectivity_reference_overlap",reference_name)
            elif scaffold_id: overlaps.append({"compound_id":record["compound_id"],"canonical_isomeric_smiles":canonical,"overlap_type":"scaffold","reference_dataset":reference_name,"reference_compound_id":scaffold_id,"exact_overlap":False,"connectivity_overlap":False,"scaffold_overlap":True})
        if exclusion_overlap: exclude(exclusion_overlap[0], exclusion_overlap[1]); continue
        flags={"pains_flag":standardized["pains_alert_count"]>0,"brenk_flag":standardized["brenk_alert_count"]>0,"extreme_molecular_weight_flag":standardized["molecular_weight"]<100 or standardized["molecular_weight"]>800,"extreme_logp_flag":standardized["logp"]<-3 or standardized["logp"]>7,"extreme_tpsa_flag":standardized["tpsa"]>250,"excessive_rotatable_bonds_flag":standardized["rotatable_bonds"]>20,"unusual_formal_charge_flag":abs(standardized["formal_charge"])>3}
        rows.append({**base,"original_smiles":original,**standardized,**flags,"structural_review_flags":";".join(k for k,v in flags.items() if v),"track1_scaffold_overlap":standardized["murcko_scaffold_smiles"] in references["track1_ki_regression"]["murcko_scaffold_smiles"],"locked_external_scaffold_overlap":standardized["murcko_scaffold_smiles"] in references["locked_external"]["murcko_scaffold_smiles"]})
    prepared=pd.DataFrame(rows); excluded=pd.DataFrame(exclusions,columns=["compound_id","original_smiles","exclusion_reason","exclusion_detail","supplier_id","catalogue_id"]); overlap_frame=pd.DataFrame(overlaps,columns=["compound_id","canonical_isomeric_smiles","overlap_type","reference_dataset","reference_compound_id","exact_overlap","connectivity_overlap","scaffold_overlap"])
    descriptor_columns=["molecular_weight","exact_molecular_weight","logp","tpsa","hbd","hba","rotatable_bonds","ring_count","aromatic_ring_count","fraction_csp3","formal_charge","qed"]
    descriptor_summary=prepared[descriptor_columns].describe().T.reset_index(names="descriptor") if not prepared.empty else pd.DataFrame(columns=["descriptor","count","mean","std","min","25%","50%","75%","max"])
    atomic_csv(prepared,output_dir/"prepared_library.csv");atomic_csv(excluded,output_dir/"excluded_compounds.csv");atomic_csv(overlap_frame,output_dir/"structure_overlaps.csv");atomic_csv(descriptor_summary,output_dir/"descriptor_summary.csv")
    alert_columns=["pains_flag","brenk_flag","extreme_molecular_weight_flag","extreme_logp_flag","extreme_tpsa_flag","excessive_rotatable_bonds_flag","unusual_formal_charge_flag"];counts=Counter(excluded.exclusion_reason) if not excluded.empty else Counter(); overlap_counts=lambda reference,kind:int(((overlap_frame.reference_dataset==reference)&(overlap_frame.overlap_type==kind)).sum()) if not overlap_frame.empty else 0
    data_outputs=[output_dir/name for name in ("candidate_pool.csv","prepared_library.csv","excluded_compounds.csv","structure_overlaps.csv","descriptor_summary.csv")];output_checksums={path.name:sha256(path) for path in data_outputs};provenance_ok=not source[list({"supplier_id","supplier_name","catalogue_id","availability","library_name","library_release","source_file"})].eq("").any().any();retained_overlap=False
    if len(prepared):
        retained_overlap=any(prepared.canonical_isomeric_smiles.isin(maps["canonical_isomeric_smiles"]).any() or prepared.inchikey.isin(maps["inchikey"]).any() or prepared.inchikey_connectivity.isin(maps["inchikey_connectivity"]).any() for maps in references.values())
    screening_ready=bool(len(prepared)>0 and provenance_ok and not retained_overlap and not prepared[["canonical_isomeric_smiles","inchikey","inchikey_connectivity"]].isna().any().any() and all(sha256(path)==output_checksums[path.name] for path in data_outputs))
    manifest={"schema_version":1,"study_id":output_dir.name,"source_library":source.library_name.iloc[0] if len(source) else None,"source_file":source.source_file.iloc[0] if len(source) else None,"supplier":source.supplier_name.iloc[0] if len(source) else None,"library_release":source.library_release.iloc[0] if len(source) else None,"availability_interpretation":"Listed in the downloaded in-stock snapshot; current stock is not guaranteed.","input_count":len(source),"invalid_structure_count":counts["invalid_or_empty_smiles"]+counts["unresolved_mixture"],"duplicate_id_count":counts["duplicate_supplier_id"],"duplicate_structure_count":counts["duplicate_standardized_structure"],"track1_exact_overlap_count":overlap_counts("track1_ki_regression","exact"),"track1_connectivity_overlap_count":overlap_counts("track1_ki_regression","connectivity"),"locked_external_exact_overlap_count":overlap_counts("locked_external","exact"),"locked_external_connectivity_overlap_count":overlap_counts("locked_external","connectivity"),"track1_scaffold_overlap_count":int(prepared.track1_scaffold_overlap.sum()) if len(prepared) else 0,"locked_external_scaffold_overlap_count":int(prepared.locked_external_scaffold_overlap.sum()) if len(prepared) else 0,"structural_alert_counts":{column:int(prepared[column].sum()) if len(prepared) else 0 for column in alert_columns},"final_usable_count":len(prepared),"input_checksum":sha256(input_path),"output_checksums":output_checksums,"git_commit":subprocess.run(["git","rev-parse","HEAD"],cwd=root,text=True,capture_output=True,check=True).stdout.strip(),"generated_timestamp":datetime.now(timezone.utc).isoformat(),"warnings":["Supplier availability means listed in the downloaded snapshot and does not guarantee current stock.","Structural alerts are review cautions, not toxicity conclusions.","No ML inference, docking, or shortlist generation was performed."],"screening_ready":screening_ready}
    schema={"candidate_pool_columns":list(source.columns),"prepared_library_columns":list(prepared.columns),"excluded_columns":list(excluded.columns),"overlap_columns":list(overlap_frame.columns)};atomic_json(schema,output_dir/"preparation_schema.json");atomic_json(manifest,output_dir/"preparation_manifest.json")
    summary=f"# Prospective library preparation\n\n- Input: {len(source):,}\n- Usable: {len(prepared):,}\n- Screening ready: {screening_ready}\n- ML inference: not run\n- Docking: not run\n- Shortlist: not created\n\nAvailability means listed in the downloaded snapshot; current stock is not guaranteed.\n";atomic_text(summary,output_dir/"preparation_summary.md")
    checksum_outputs=[*data_outputs,output_dir/"preparation_manifest.json",output_dir/"preparation_schema.json",output_dir/"preparation_summary.md"]
    atomic_json({path.name:sha256(path) for path in checksum_outputs},output_dir/"preparation_checksums.json")
    if strict and not screening_ready: raise RuntimeError("Strict preparation failed; inspect exclusion and readiness reports")
    return manifest
