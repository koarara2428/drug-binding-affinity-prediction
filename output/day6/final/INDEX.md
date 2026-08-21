# 최종 산출물 안내

`output/day6/` 전체(35개 파일)에서 결론을 이해하고 모델을 쓰는 데 필요한 핵심만 12개로 추렸다. 나머지 원본 파일은 `output/day6/figures/`, `output/day6/models/`에 그대로 남아있다(삭제 안 함).

## 실제로 쓸 것

| 파일 | 내용 |
|---|---|
| `final_model.pkl` | **최종 모델.** Ridge(9피처, StandardScaler), train+test 전체(1,825행)로 재학습. `joblib.load()`로 불러오면 `{"model", "feature_columns"}` 딕셔너리 |
| `final_model_coefficients.csv` | 최종 모델의 표준화 회귀계수. `logp_pi_interaction`이 압도적 1위 |

## 왜 이 모델을 골랐는지 (근거)

| 파일 | 내용 | 관련 README |
|---|---|---|
| `model_comparison.csv` / `.png` | Linear Regression / RandomForest / XGBoost 3모델 1차 비교(단일 분할) — `.png`는 README에 실제 삽입된 그림 | 모델 비교 |
| `cross_validation_results.csv` | 위 3모델의 5-fold 교차검증 — 성능 차이가 통계적으로 무의미했다는 근거 | 모델 비교 |
| `next_steps_summary.csv` | 상호작용항 / `protein_id` 그룹분할 / Ridge(7·9피처) 실험 결과(단일 분할) | 실험 기록 |
| `next_steps_cv_results.csv` | 위 실험들의 교차검증 — **Ridge 9피처가 유의미하게 낫다**는 최종 근거(R² 0.5714±0.0681) | 모델 비교, 실험 기록 |

## 모델의 한계 (같이 봐야 할 것 — 가장 중요)

| 파일 | 내용 | 관련 README |
|---|---|---|
| `residual_distribution.png` | (baseline 7피처) 오차 분포 — 대부분 작지만 꼬리가 김 | 핵심 발견 |
| `predicted_vs_residual.png` | (baseline 7피처) 예측값별 오차 패턴 | 핵심 발견 |
| `linear_regression_top10_errors.csv` | (baseline 7피처) 오차 최대 10건 — 회귀-평균 편향의 증거 | 핵심 발견 |
| `final_model_residual_distribution.png` | **(최종 모델) 오차 분포** | 핵심 한계 |
| `final_model_predicted_vs_residual.png` | **(최종 모델) 예측값별 오차 패턴 — README에 실제 삽입된 그림** | 핵심 한계 |
| `final_model_top10_errors.csv` | **(최종 모델) 오차 최대 10건 — R²는 올랐지만 이 편향은 거의 그대로 남음(corr 0.773→0.736)** | 핵심 한계 |

## 데이터 이해 (핵심 EDA)

| 파일 | 내용 | 관련 README |
|---|---|---|
| `correlation_heatmap_final.png` | 전처리 최종본 기준 전체 상관관계 | 핵심 발견 |
| `binding_affinity_histogram.png` | 타깃 분포 (평균 6.52, 오른쪽 치우침) | 프로젝트 요약 |

**참고**: README.md는 2026-08-21 채용 담당자용 포트폴리오 형식("Claude Code로 수행한 신약 탐색 데이터 분석")으로 전면 재구성됐다(한글). 위 표는 그 최신 구조 기준이다. 새로 추가된 `docs/CLAUDE_WORKFLOW.md`, `docs/biological_context.md`, `docs/next_experiment.md`, `output/day6/final/FINAL_RESULTS.md`도 함께 참고.

---
전체 과정과 수치 해석은 프로젝트 루트의 `README.md`를 참고. 시간순 작업 기록은 `docs/WORKFLOW.md` 참고.
