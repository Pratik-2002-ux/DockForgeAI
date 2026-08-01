def metric_columns(frame): return [c for c in ("pr_auc","roc_auc","mcc","brier_score","expected_calibration_error") if c in frame]
