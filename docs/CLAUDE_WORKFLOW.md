# Claude Code Workflow

> **Representative workflow.** This document summarizes how Claude Code was used across this project, based on the actual sequence of analysis steps and decisions recorded in [`WORKFLOW.md`](WORKFLOW.md). It is not a verbatim chat log — no such log file exists in this repository — but a structured reconstruction of the real prompt → code → validation → decision cycle that was followed at each stage.

## Objective

Show how Claude Code was used as a tool inside a researcher-directed analysis: Claude Code generated code and surfaced results; the researcher set the problem, judged data quality, evaluated statistical evidence, and decided what to do next.

## 1. Problem Definition

- **Instruction**: Define what the dataset supports (a regression problem), confirm the target (`binding_affinity`) and a plausible feature set, and state a concrete success criterion (a model usable as a pre-experiment screening aid).
- **Claude Code action**: Ran a structured feasibility check (row/column count, candidate continuous targets, missing-value ratios, target-leakage screen) before any modeling.
- **Researcher validation**: Confirmed the problem was well-posed for regression and that the stated goal (screening aid, not a production tool) was realistic for the dataset size and feature quality.
- **Decision**: Proceed with `binding_affinity` as the target; document the goal in `docs/problem_definition.md`.

## 2. Data Inspection

- **Instruction**: Inspect the raw file's structure, dtypes, and value ranges before touching anything.
- **Claude Code action**: Produced `shape`, `dtypes`, `nunique`, missing-value ratios, and a raw-string re-read to check for disguised missing tokens (`"?"`, `"N/A"`, blanks) and numbers stored as text (commas, units).
- **Researcher validation**: Reviewed the printed diagnostics directly; confirmed no disguised missing values or type issues existed, so no extra cleanup step was needed for those two categories.
- **Decision**: Move to a formal data-quality assessment rather than assuming the raw file was clean.

## 3. Data Quality Assessment

- **Instruction**: Check for duplicates, invalid values, and any column that might leak the target.
- **Claude Code action**: Checked full-row duplicates (0 found), scanned numeric columns for implausible values, and computed the relationship between `active` and `binding_affinity`.
- **Researcher validation**: Identified that `active` is an exact binarization of `binding_affinity` at a threshold of 7.0 (active=0 max 6.996, active=1 min 7.002) — a textbook case of target leakage — and that one `polar_surface_area` value (-24.65) was physically impossible (surface area cannot be negative).
- **Decision**: Drop `active` and `compound_id` (pure identifier) as features; reclassify the invalid `polar_surface_area` value as missing; drop rows with missing values in `logp`, `polar_surface_area`, or `hydrophobicity` (174 rows) plus the one invalid-value row. Dataset reduced from 2,000 to 1,825 rows, 15 columns.

## 4. Exploratory Data Analysis

- **Instruction**: Characterize the target distribution and each feature's relationship to it before choosing features.
- **Claude Code action**: Computed target summary statistics, a full correlation matrix, per-feature distribution plots, and an outlier scan (IQR-based).
- **Researcher validation**: Read the correlation ranking directly — only `logp`, `protein_pi`, and (via a derived column) `logp_pi_interaction` showed meaningful correlation with the target; the originally expected driver, `molecular_weight`, did not. Flagged two multicollinear pairs (`logp`↔`logp_pi_interaction`, `protein_length`↔`mw_ratio`) for the feature-selection stage.
- **Decision**: Select a 7-feature baseline set excluding one member of each multicollinear pair, to keep the first model's coefficients interpretable.

## 5. Modeling

- **Instruction**: Start with an interpretable baseline, then compare against more flexible models under identical train/test conditions.
- **Claude Code action**: Trained Linear Regression, RandomForest, and XGBoost on the same 8:2 split (`random_state=42`) and reported R²/MAE/RMSE for each.
- **Researcher validation**: On the single split, RandomForest and XGBoost appeared to outperform Linear Regression — but this observation was treated as provisional, not conclusive.
- **Decision**: Do not select a "best" model from a single split; proceed to formal validation before drawing conclusions.

## 6. Model Validation

- **Instruction**: Re-test the single-split ranking with 5-fold cross-validation before trusting it.
- **Claude Code action**: Ran `KFold(5, shuffle=True, random_state=42)` cross-validation on all three models and reported R² mean ± std per model.
- **Researcher validation**: Applied a simple decision rule (a performance gap smaller than the larger model's standard deviation is not meaningful) and found **all three models were statistically indistinguishable** — the single-split "advantage" for RandomForest/XGBoost fell inside normal fold-to-fold variability.
- **Decision**: Since raw model choice did not move performance, investigate feature representation instead — specifically, whether Ridge regularization could safely reintroduce the two features dropped for multicollinearity.

## 7. Residual Analysis

- **Instruction**: Do not stop at an aggregate R²; check where the model is wrong, not just how wrong on average.
- **Claude Code action**: Computed residuals on the held-out test set, plotted their distribution and their relationship to predicted values, listed the 10 largest-error rows, and correlated residuals against the actual target.
- **Researcher validation**: All 10 largest errors occurred at the extremes of the target range (very high or very low actual affinity); `corr(residual, actual)` was 0.77 for the baseline model — a strong regression-to-the-mean bias. After the Ridge (9-feature) model raised overall R², the same check was repeated on that model: the correlation only dropped to 0.74, meaning the bias persisted.
- **Decision**: Document this as a critical, unresolved limitation rather than treating the R² improvement as a full fix.

## 8. Researcher Decision Points

These decisions were made by the researcher based on Claude Code's output, not automated by it:

- Which columns constitute target leakage vs. legitimate features
- How to treat missing and invalid values (drop vs. impute), and at what threshold
- Which features to include in the interpretable baseline vs. the regularized final model
- Whether a single-split performance difference was worth acting on (it was re-tested via cross-validation first)
- Whether the final model's residual bias was acceptable for the stated use case (it was not, for final candidate ranking)
- How to scope the model's practical use (preliminary filter, not final ranking tool)

## 9. Lessons Learned

- A single train/test split is not sufficient evidence to rank models; cross-validation reversed the initial ranking conclusion.
- Removing correlated features to keep a model interpretable is a valid but not universal choice — regularization (Ridge) let the same information back in without destabilizing coefficients, and materially improved performance (CV R² 0.455 → 0.571).
- An improved aggregate metric does not guarantee a specific weakness is fixed. The high-affinity underprediction bias was checked explicitly on the final model, not assumed to be resolved.
