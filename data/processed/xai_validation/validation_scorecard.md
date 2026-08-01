# DockForgeAI XAI Validation Scorecard

Qualification: **benchmark-supported with cautions**

No overall percentage is reported.

## Dataset Integrity

- Status: **pass**
- Evidence: Locked manifest and checksums
- Limitation: Internal checks do not prove labels

## Duplicate Control

- Status: **pass**
- Evidence: Curated duplicate exclusion report
- Limitation: Source records may share analog series

## Scaffold Leakage

- Status: **pass**
- Evidence: Scaffold-safe split
- Limitation: Scaffold separation is not external validation

## Source Bias

- Status: **fail**
- Evidence: Source-only classifier
- Limitation: Source and label provenance are coupled

## Adversarial Separability

- Status: **caution**
- Evidence: Chemical-domain classifier
- Limitation: In-sample adversarial estimate

## Permutation Significance

- Status: **pass**
- Evidence: 100 permutations
- Limitation: Significance does not remove source bias

## Dummy Baseline Superiority

- Status: **pass**
- Evidence: Locked-test dummy baselines
- Limitation: Same benchmark limitations apply

## Nearest Neighbour Superiority

- Status: **caution**
- Evidence: Tanimoto and descriptor neighbours
- Limitation: Analog memorization can inflate results

## Matched Subset Robustness

- Status: **not_evaluated**
- Evidence: No-replacement descriptor matching
- Limitation: No pairs met the configured caliper

## Decoy Subset Robustness

- Status: **pass**
- Evidence: Cached decoy sensitivity
- Limitation: DUD-E decoys remain computational

## Calibration Quality

- Status: **caution**
- Evidence: Locked-test and domain calibration
- Limitation: Sparse out-of-domain subgroup

## Applicability Domain Reliability

- Status: **caution**
- Evidence: 263 in-domain / 37 out-of-domain
- Limitation: AD is similarity-based

## Seed Stability

- Status: **pass**
- Evidence: 5 seeds
- Limitation: Same fixed split

## Explanation Sanity

- Status: **caution**
- Evidence: Randomization/removal/agreement controls
- Limitation: Non-causal explanations

## External Validation Status

- Status: **not_evaluated**
- Evidence: No user-supplied independent dataset
- Limitation: No external evidence
