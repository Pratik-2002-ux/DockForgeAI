# XAI features

`build_xai_features.py` creates 33 finite curated descriptors, 2,048 radius-2 Morgan bits, 167 MACCS keys, and a 2,248-column fused sparse matrix. Identifiers, labels, sources, scaffolds, activity values, splits, and provenance are prohibited model inputs. Row order is fixed by `feature_index.csv` and tied to the dataset checksum.

```bash
conda run -n dockforge python scripts/build_xai_features.py --dataset data/processed/xai_dataset/hivpr_xai_2000_compounds.csv --overwrite
```
