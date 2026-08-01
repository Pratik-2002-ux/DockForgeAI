# XAI limitations

- DUD-E decoys are computational controls, not experimentally proven inactive molecules.
- ChEMBL/DUD-E source differences can create shortcut signals and inflate apparent classification.
- Good classification does not prove binding, inhibition, therapeutic value, or clinical relevance.
- SHAP/native importance explains model behaviour, not biological causality.
- Fingerprint hashing causes collisions.
- Applicability-domain status measures training support, not truth.
- No molecular generation, optimization, docking, active-learning proposal, wet-lab planning, or experimental validation occurs in this milestone.
# Additional validation limitations

ChEMBL actives and DUD-E computational decoys can be separable by source-specific chemistry. DUD-E decoys are not experimentally verified inactive compounds. High benchmark scores do not prove binding, therapeutic activity, or generalization to independent chemistry. Independent external validation has not been completed and experimental validation is outside scope.
