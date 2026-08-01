# XAI permutation tests

Training labels are permuted while locked validation and test labels remain unchanged for scoring. The empirical p-value is `(scores >= observed + 1) / (repetitions + 1)`. The report always records the repetition count. Statistical separation from permuted labels does not eliminate source bias or establish therapeutic activity.
