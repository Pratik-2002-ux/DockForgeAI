# XAI validation

The offline validation suite tests whether the locked benchmark is robust to source shortcuts, label permutations, trivial baselines, chemical similarity, scaffold frequency, physicochemical matching, decoy selection, calibration, seed variation, and explanation controls. It consumes only repository-cached data and never docks, generates molecules, or tunes the locked model on test data.

Software tests establish deterministic implementation behaviour. They do not establish biological truth. The current qualification is **benchmark-supported with cautions**, not externally or experimentally validated.

Run `conda run --no-capture-output -n dockforge python scripts/run_xai_validation_suite.py --all --permutations 100 --bootstrap-iterations 500 --decoy-seeds 2026 2027 2028 2029 2030 2031 2032 2033 2034 2035 --offline --overwrite`.
