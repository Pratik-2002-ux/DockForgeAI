# XAI modelling

Logistic regression, random forest, extra trees, and histogram gradient boosting are evaluated on descriptor, Morgan, and fused features across seeds 2026–2030. Existing Bemis–Murcko splits are never reshuffled: training fits models, validation selects/calibrates, and test remains locked. Validation PR-AUC is primary.

High benchmark scores may reflect the ChEMBL-active versus DUD-E-decoy source distinction and must not be interpreted as binding or therapeutic performance.
