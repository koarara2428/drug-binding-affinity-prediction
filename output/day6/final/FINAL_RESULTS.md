# Final Results

## Dataset

화합물-단백질 조합 2,000건, 17컬럼(원본). 품질 검토 후: 1,825행, 15컬럼.

## Target

`binding_affinity` (연속형)

## Best Model

Ridge 회귀, 9피처, `StandardScaler` + `alpha=1.0`

## Performance

| 지표 | 값 |
|---|---|
| CV R² (5-fold, mean ± std) | 0.5714 ± 0.0681 |
| Test R² (단일 8:2 분할) | 0.5533 |
| Test MAE | 0.3700 |
| Test RMSE | 0.8759 |

비교 기준(baseline) — Linear Regression, 7피처: CV R² = 0.4552 ± 0.0577.

## Key Finding

`logp`와 `protein_pi`의 상호작용만이 `binding_affinity`와 개별적으로 의미 있는 상관관계를 보인 피처였습니다. 다중공선성 때문에 처음에 제외했던 두 피처(`logp_pi_interaction`, `mw_ratio`)를 수동 제외 대신 Ridge 정규화로 다시 포함시키자 교차검증 R²가 0.455에서 0.571로 상승했고, 이는 통계적으로 유의미한 개선이었습니다(차이가 두 모델의 fold 간 표준편차보다 큼).

## Model Limitation

baseline과 최종 모델 모두 고결합력 화합물을 체계적으로 과소예측합니다. held-out test set 기준 `corr(residual, 실제값)`은 baseline이 0.77, 최종 모델이 0.74였습니다 — 전체 R² 개선이 이 편향을 해소하지 못했습니다. 최종 모델에서 실제 결합력 상위 10% 화합물은 평균 +0.96(residual)만큼 과소예측됐고, 중간 80% 구간은 편향이 거의 없었습니다(+0.01).

## Practical Interpretation

이 모델은 **최종 후보 순위화 모델**보다는 **1차 필터링 도구**에 더 적합합니다 — 모델이 약한 결합체일 것으로 확신하는 후보를 우선순위에서 제외하는 용도입니다. 결합력이 가장 높은 화합물에서 특히 부정확하기 때문에, 예측 순위를 기준으로 최종 후보를 선택하면 실제로 강력한 후보를 놓칠 위험이 있습니다.

## Next Experiment

제안하는 후속 작업(고결합력 가중 회귀, 분위수 회귀, XGBoost 튜닝, 구조 기반 피처, 외부 검증)은 [`docs/next_experiment.md`](../../../docs/next_experiment.md)를 참고하세요.
