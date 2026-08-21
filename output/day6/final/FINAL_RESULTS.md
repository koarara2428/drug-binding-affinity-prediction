# Final Results

## Dataset

2,000 compound-protein pairs, 17 columns (raw). After quality filtering: 1,825 rows, 15 columns.

## Target

`binding_affinity` (continuous)

## Best Model

Ridge Regression, 9 features, `StandardScaler` + `alpha=1.0`

## Performance

| Metric | Value |
|---|---|
| CV R² (5-fold, mean ± std) | 0.5714 ± 0.0681 |
| Test R² (single 8:2 split) | 0.5533 |
| Test MAE | 0.3700 |
| Test RMSE | 0.8759 |

Baseline for comparison — Linear Regression, 7 features: CV R² = 0.4552 ± 0.0577.

## Key Finding

`logp` and its interaction with `protein_pi` are the only features with a meaningful individual correlation to `binding_affinity`. Reintroducing two features (`logp_pi_interaction`, `mw_ratio`) that were initially excluded for multicollinearity — using Ridge regularization instead of manual exclusion — raised cross-validated R² from 0.455 to 0.571, a statistically meaningful improvement (gap larger than either model's fold-to-fold standard deviation).

## Model Limitation

Both the baseline and final models systematically underpredict high-affinity compounds. On the held-out test set, `corr(residual, actual value)` was 0.77 for the baseline and 0.74 for the final model — the overall R² improvement did not resolve this bias. In the final model, the top 10% of compounds by actual affinity were underpredicted by an average of +0.96 (residual), while the middle 80% showed near-zero bias (+0.01).

## Practical Interpretation

The model is better positioned as a **preliminary filtering tool** — for deprioritizing candidates it confidently predicts as weak binders — rather than a **final candidate-ranking model**. Because it specifically struggles on the highest-affinity compounds, using its predicted ranking to select a final shortlist risks excluding genuinely strong candidates.

## Next Experiment

See [`docs/next_experiment.md`](../../../docs/next_experiment.md) for proposed follow-ups: high-affinity weighted regression, quantile regression, XGBoost tuning, structure-derived features, and external validation.
