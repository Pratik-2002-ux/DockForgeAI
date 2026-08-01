# Reporting

Generate both formats:

```bash
conda run -n dockforge python scripts/generate_final_report.py --overwrite
```

Use `--pdf-only` or `--html-only`, `--strict` to reject missing optional inputs, and `--output-prefix NAME` to change filenames. Existing outputs are protected unless `--overwrite` is supplied. `--open` is the only option that launches a macOS viewer.

The orchestrator aggregates JSON/CSV outputs, renders deterministic RDKit molecule images, generates matplotlib figures, creates PDF and offline HTML reports, validates non-empty files, and writes `data/processed/final_reports/report_manifest.json`. Full failures are logged under `logs/reporting/`; optional molecule/figure failures become manifest warnings.

The reports label Vina scores as approximate, the demo screen as software validation, descriptors as distinct from full ADMET, and predictions as requiring experiments.
