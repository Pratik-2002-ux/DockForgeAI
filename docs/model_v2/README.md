# DockForgeAI v2

DockForgeAI v2 is an isolated, experimental modelling framework. DockForgeAI
v1 remains the validated platform and is not changed by v2 training code.

## Current status

Readiness only. The repository contains cached ChEMBL activity records, but it
does not contain an approved independent BindingDB training export. The local
PubChem records are the locked external benchmark and must never be used for
training, calibration, threshold selection, or model selection.

Run the offline readiness audit:

```bash
conda run -n dockforge-v2 python scripts/model_v2/audit_v2_readiness.py
```

Training must remain disabled while the generated report has
`training_allowed: false`.
