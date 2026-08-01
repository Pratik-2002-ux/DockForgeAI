# Developer guide

Keep scientific modules typed, documented, and independent of Streamlit. Resolve paths through `pathlib.Path` and `config/settings.yaml`; do not write private absolute paths into reports. Use custom exceptions at API boundaries and logging for report failures. Never run docking in unit tests.

Validation:

```bash
conda run -n dockforge python -m compileall app src scripts tests
conda run -n dockforge pytest --collect-only -q
conda run -n dockforge pytest -q
conda run -n dockforge python scripts/generate_final_report.py --overwrite
```

Preserve existing outputs and the corrected 0.518 Å RMSD baseline. Add focused tests for new behavior, do not weaken assertions, and keep generated caches/logs out of Git. External ADMET integrations must retain source, version or assay, units, and retrieval date. Molecular dynamics and experimental validation are outside current scope.
