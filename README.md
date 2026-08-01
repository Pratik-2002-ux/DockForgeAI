# DockForgeAI

DockForgeAI is a reproducible, resource-aware molecular docking, virtual-screening, interaction-analysis, and hit-prioritization workflow with an interactive Streamlit dashboard and portable scientific reports.

![Python 3.12](https://img.shields.io/badge/Python-3.12-blue) ![Tests](https://img.shields.io/badge/tests-57%20passing-brightgreen) ![License](https://img.shields.io/badge/license-not%20yet%20specified-lightgrey)

## Scientific overview

The project turns prepared protein and ligand structures into validated redocking poses, geometric interaction summaries, small-library screening rankings, RDKit descriptor profiles, and consensus prioritization. It exists to keep those steps auditable and practical on an Apple Silicon laptop with 8 GB RAM.

DockForgeAI does not turn docking scores into experimental claims. Vina scores are approximate ranking values, RDKit descriptors are not full ADMET predictions, and every predicted hit requires experimental validation.

## Key features

- Protein retrieval, inspection, cleaning, and receptor preparation
- Ligand chemistry correction, conformer handling, and PDBQT preparation
- AutoDock Vina redocking, corrected multi-pose SDF export, and crystallographic RMSD validation
- Protein–ligand interaction comparison and catalytic Asp25 analysis
- Resumable, low-memory virtual screening and consensus hit prioritization
- Interactive py3Dmol pose viewer and RDKit 2D molecule rendering
- Offline PDF/HTML reports, manifest, download bundles, Docker packaging, and CI

## Architecture and workflow

```text
Input structures → protein/ligand preparation → docking box → Vina redocking
       → RMSD validation → interaction analysis → demonstration screening
       → descriptors/alerts → consensus prioritization → PDF/HTML reports
```

The `src/` tree contains testable domain modules, `scripts/` contains command-line entry points, `app/` provides Streamlit, `config/settings.yaml` centralizes paths and parameters, and `data/processed/` retains validated outputs. See [architecture](docs/ARCHITECTURE.md).

## Validated 1HSG/MK1 case study

The HIV-1 protease reference case retains chains A/B and the catalytic Asp25 dyad. The corrected redocking result is:

- Top-ranked Vina score: **−11.575 kcal/mol**
- Crystallographic heavy-atom RMSD: **0.518 Å**
- Rank 1 passes the 2 Å criterion
- Interpretation: **Validated strongly**

The three-compound virtual screen completed 3/3 docking runs; compound 3 ranked first at −2.848 kcal/mol. This validates software execution only. It is not a biological discovery and does not establish any compound as an inhibitor.

## Screenshots

Add project-owned screenshots at `docs/images/dashboard.png`, `docs/images/viewer.png`, and `docs/images/export-center.png`. These placeholders intentionally avoid external assets.

## Installation with Conda

```bash
conda env create -f environment.yml
conda activate dockforge
```

For an existing environment, use `conda env update -n dockforge -f environment.yml --prune`. Detailed platform notes are in [installation](docs/INSTALLATION.md).

## Quick start

Launch the dashboard:

```bash
conda run -n dockforge streamlit run app/streamlit_app.py
```

Generate final reports from existing outputs:

```bash
conda run -n dockforge python scripts/generate_final_report.py --overwrite
```

Open reports through the **Reports & Exports** page or under `data/processed/final_reports/`.

## Docker

```bash
docker compose up --build
```

Then open <http://localhost:8501>. Docker uses a Linux micromamba image and persistent `data`, `reports`, and `logs` mounts. Apple Silicon package availability and validation status are documented in [Docker notes](docs/DOCKER.md).

## Inputs and outputs

Supported inputs include PDB receptors, SDF/SMILES ligands, PDBQT docking inputs/poses, and CSV libraries or external ADMET tables with declared provenance. Generated outputs include PDB/PDBQT/SDF structures, CSV/JSON/Markdown analyses, PNG/SVG figures, PDF/HTML reports, and a JSON report manifest.

```text
data/processed/
├── docking_results, docking_tables, docking_poses
├── interaction_reports
├── screening_tables, screening_reports
├── prioritization_reports
├── report_assets/{figures,molecules}
└── final_reports/{pdf,html,report_manifest.json}
```

## Testing and reproducibility

```bash
conda run -n dockforge python -m compileall app src scripts tests
conda run -n dockforge pytest -q
```

Seeds, Vina settings, directory paths, component weights, and catalytic residues live in `config/settings.yaml`. Reports retain project-relative source provenance and the current Git hash when available. Expensive docking is excluded from report and CI tests.

## Hardware constraints

The workflow is intentionally serialized and chunked for an 8 GB Apple Silicon system. It has no GPU requirement. Larger libraries need deliberate batching and storage planning.

## Scientific limitations

- Docking scores are approximate; they are not experimental binding free energies.
- Geometric interactions depend on protonation and distance/angle assumptions.
- Descriptor and alert calculations do not replace measured ADMET.
- The demo screen is too small for discovery claims.
- Molecular dynamics, free-energy calculations, biochemical assays, and cellular validation are outside the current scope.

## Roadmap

Potential next work includes larger benchmark datasets, richer provenance schemas, calibrated uncertainty, and experimentally anchored validation. No roadmap item should weaken the current reproducible baseline.

## Citation

A formal citation and archival DOI have not yet been assigned. Until then, cite the repository name, version, commit hash, and access date.

## License

No license has yet been specified. Do not assume redistribution or commercial-use rights until a license file is added.

## Contributing

Read the [developer guide](docs/DEVELOPER_GUIDE.md), keep changes scoped, add tests, preserve the 0.518 Å reference result, and run the full compile/test commands before proposing changes.

Further guides: [user guide](docs/USER_GUIDE.md) · [reporting](docs/REPORTING.md) · [Docker](docs/DOCKER.md)
