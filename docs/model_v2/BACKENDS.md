# Optional v2 model backends

Random Forest uses the existing scikit-learn environment. XGBoost, LightGBM, CatBoost, Chemprop, PyTorch Geometric graph models, ONNX, and TorchScript are explicit optional backends. Missing backends are reported as `not_evaluated`; v2 never substitutes a different algorithm under the requested name.

Current frozen environment availability is recorded by the experiment manifest. Install optional backends in a separate v2 environment rather than altering the frozen v1 environment.
