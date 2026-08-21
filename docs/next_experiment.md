# Next Experiment

This document lists proposed follow-up analyses based on the limitations identified in this project. **None of the items below have been run.** They are proposals for future work, not completed experiments.

## Hypothesis

The final model (Ridge, 9 features) systematically underpredicts high-affinity compounds (`corr(residual, actual)` = 0.74 on the held-out test set, barely improved from the 7-feature baseline's 0.77). Two candidate explanations were not yet distinguished:

1. **Loss function mismatch**: standard squared-error loss optimizes for average accuracy, which favors the dense middle of the target distribution over the sparse high-affinity tail.
2. **Insufficient feature information**: the available molecular/protein descriptors may not capture what actually drives unusually strong binding, regardless of model choice (RandomForest and XGBoost showed the same underlying bias in earlier checks).

Both may be true simultaneously; the proposed experiments below are designed to tell them apart.

## Proposed Experiments

- **High-affinity weighted regression** — reweight training samples (e.g., inversely to local target density, or proportionally to `|y - mean|`) so errors on high-affinity compounds cost more during training. Tests explanation (1) directly.
- **Quantile regression** — train a model to predict an upper quantile (e.g., 90th percentile) of binding affinity instead of the mean, which does not have the same incentive to shrink toward the average.
- **XGBoost hyperparameter tuning** — the XGBoost model used in this project used untuned, conservative defaults (`max_depth=3, learning_rate=0.05`). A tuned search (e.g., `GridSearchCV`) was not attempted and might close some of the gap, though it is unlikely to fix a bias already observed across three different model families.
- **Structure-derived features** — descriptors beyond basic molecular/protein properties (e.g., docking scores, 3D interaction features) would test explanation (2): if the bias persists even with richer features, the limitation is more likely in the loss function than the feature set.
- **External validation** — evaluate the current model (or any improved version) against a dataset not used in this project, ideally one with experimentally measured affinities rather than simulated ones.

## Evaluation Metrics

Beyond the metrics already used in this project (R², MAE, RMSE), the following would specifically target the high-affinity underprediction issue:

- **Spearman rank correlation** — measures whether relative ranking is preserved even if absolute values are biased, which matters more than raw R² for a prioritization use case.
- **Top-K recall / precision** — for a fixed shortlist size K, what fraction of the true top-K highest-affinity compounds does the model's top-K predicted list actually contain? This is arguably the most direct metric for a screening tool, since the practical question is "does it find the good candidates," not "is its average error small."

Regression accuracy alone (R²/MAE/RMSE) does not answer whether a screening model retrieves the compounds that matter most — that requires ranking-focused metrics like the two above.
