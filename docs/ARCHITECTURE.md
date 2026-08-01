# Architecture

`src/proteins` and `src/ligands` prepare inputs; `src/docking` generates boxes, runs/parses Vina, exports poses, and validates RMSD; `src/analysis` compares interactions; `src/screening` manages resumable batches; `src/admet` and `src/prioritization` calculate and rank features; `src/visualization` renders structures; and `src/reports` aggregates immutable outputs into figures and portable reports.

The Streamlit layer in `app/` calls read-only loaders and bounded report helpers. Configuration is centralized in `config/settings.yaml`. Reports preserve source provenance as project-relative paths. CLI scripts orchestrate modules but do not duplicate scientific logic.

The design favors one-process, low-memory execution on 8 GB systems. CI compiles and tests small fixtures without docking or internet access.
