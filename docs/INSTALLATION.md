# Installation

## Conda

DockForgeAI targets Python 3.12 and uses conda-forge/bioconda for RDKit, Vina, Gemmi, and Open Babel.

```bash
conda env create -f environment.yml
conda activate dockforge
conda run -n dockforge pytest -q
```

Update an existing environment with `conda env update -n dockforge -f environment.yml --prune`. The environment contains no GPU runtime. Confirm tools with `vina --version` and `python -c "import rdkit, meeko, gemmi, streamlit, py3Dmol"`.

macOS Apple Silicon is the validated host. Linux is targeted through the same environment definition, but solver availability can change; review any solve before accepting package substitutions.
