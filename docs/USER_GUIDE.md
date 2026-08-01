# User guide

Launch the application with:

```bash
conda run -n dockforge streamlit run app/streamlit_app.py
```

Use **Redocking** for the corrected 1HSG/MK1 validation, **3D viewer** to rotate/zoom the receptor and compare three poses, **Interactions** for contact tables, **Virtual screening** and **Hit prioritization** for the demonstration rankings, and **Reports & Exports** for portable files and ZIP bundles.

Generate reports without the UI using `conda run -n dockforge python scripts/generate_final_report.py --overwrite`. This reads existing results and does not rerun Vina. External ADMET may only be interpreted when its source, units, model/assay, and collection date are recorded. Missing provenance is reported as unavailable.

The demo screen demonstrates software operation, not biological discovery. Molecular dynamics and experimental validation are outside the application’s current scope.
