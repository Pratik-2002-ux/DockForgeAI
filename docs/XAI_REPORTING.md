# XAI reporting

```bash
conda run -n dockforge python scripts/generate_xai_report.py
```

Outputs include high-resolution PNG/SVG figures, JSON/Markdown model cards, and portable PDF/HTML reports with checksummed manifest entries. Reports state the decoy, source-shortcut, hashing, causality, applicability-domain, and experimental-validation limitations.
