# HIV-1 protease XAI benchmark dataset

Phase XAI-1 constructs and validates a benchmark of up to 2,000 unique structures. It does not dock compounds or train a model.

## Sources and labels

- **ChEMBL:** official activity API, target `CHEMBL243` (human immunodeficiency virus type 1 protease). Exact numeric IC50, Ki, Kd, or EC50 records are retained when target confidence is at least 8 where reported. `active` requires pChEMBL ≥ 6.0 or standardized potency ≤1,000 nM.
- **DUD-E:** official `HIVPR` target files. Selected `decoys_final.ism` structures receive `binary_label=inactive_or_decoy` and `label_source=dude_decoy`.

DUD-E decoys are property-matched computational controls. They are not experimentally proven inactive and must not be merged conceptually with measured ChEMBL inactivity. Measurements with pChEMBL between 5.0 and 6.0 are ambiguous and excluded from binary classification.

## Standardization and conflicts

RDKit sanitizes structures, retains the largest organic fragment, removes simple disconnected counterions, preserves stereochemistry, and produces canonical isomeric SMILES and InChIKeys. Duplicate audits use canonical SMILES, full InChIKey, and the connectivity block. Replicate ChEMBL records retain median/minimum/maximum pChEMBL and spread. Compounds are excluded when records imply both active and inactive labels or span more than 2 log units.

Decoys are deterministically selected with seed 2026 to approximate active molecular weight, LogP, HBD, HBA, rotatable bonds, TPSA, and formal charge. PAINS and Brenk alerts are recorded, not silently filtered.

## Leakage control

Bemis–Murcko scaffold groups are assigned intact to training (70%), validation (15%), or test (15%). The validator rejects any scaffold appearing in more than one split. An ordinary random split is not the primary dataset split.

## Build

```bash
conda run -n dockforge python scripts/build_hivpr_xai_dataset.py --target-count 2000 --active-count 1000 --decoy-count 1000 --resume
```

Use `--offline` after raw responses and target files are cached, `--refresh` to deliberately retrieve new source data, `--strict` for downstream automation, and `--seed` to control deterministic matching/splitting. API pages are cached after every response and downloads resume from existing validated files.

Validate independently:

```bash
conda run -n dockforge python scripts/validate_hivpr_xai_dataset.py --dataset data/processed/xai_dataset/hivpr_xai_2000_compounds.csv
```

Raw page caches, archives, partial files, and logs are excluded from Git. Final manifests retain source endpoints, release information when reported, filenames, timestamps, checksums, and project-relative references.

## Scientific and ethical limitations

Dataset composition reflects historical assay coverage and source-database biases. Numeric activities across assay contexts are not perfectly interchangeable. Decoys can introduce shortcut-learning artifacts, so future model evaluation must compare scaffold-aware and source-aware behavior. No model performance is claimed. Docking, feature engineering, machine learning, SHAP, active learning, molecular dynamics, and experimental validation are outside Phase XAI-1.
