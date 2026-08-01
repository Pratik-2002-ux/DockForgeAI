# XAI dashboard

The dashboard adds Overview, Model Performance, Global Explanations, Compound Explanation, Applicability Domain, Benchmarking, and XAI Reports pages. Pages load cached artifacts and never train automatically.

```bash
conda run -n dockforge streamlit run app/streamlit_app.py
```
# Validation page

The XAI Validation page displays cached source-bias, permutation, baseline, similarity, scaffold, matching, decoy, calibration, explanation-sanity, external-interface, scorecard, and report artifacts. It never starts training automatically and uses caution—not validated—indicators where evidence is unresolved.
