# Binding Affinity 예측 — 분석 보고서

## Executive Summary

- **목적**: 화합물-단백질 결합 친화도(`binding_affinity`)를 물리화학적 특성으로 예측해, 신약 후보를 실험 전 컴퓨터로 1차 스크리닝하는 데 참고할 모델을 만든다.
- **최종 권장 모델**: **Ridge 9피처(StandardScaler, alpha=1.0)**. 교차검증 기준 **R² = 0.5714 ± 0.0681**.
- **베이스라인**: 7피처 Linear Regression(교차검증 R² = 0.4552 ± 0.0577)을 기준선으로 별도 분리해 다뤘으며, 최종 모델은 이보다 통계적으로 유의미하게 우수하다.
- **핵심 한계**: 최종 모델도 결합력이 특히 높은(스크리닝에서 가장 중요한) 후보를 체계적으로 과소예측하는 편향이 남아 있다. 전체 설명력은 개선됐지만 이 문제는 거의 해결되지 않았다.

---

## 1. 문제 정의

`binding_affinity`를 `molecular_weight`, `logp`, `protein_pi` 등 화합물·단백질 특성으로 예측하는 회귀 문제다. 대상 사용자는 신약 후보를 스크리닝하는 연구원이며, 성공 기준은 실험 전 1차 필터링에 참고할 수 있는 수준의 예측력이다.

## 2. 데이터셋 개요

- 원본: `data/raw/drug_discovery_virtual_screening.csv`, 2,000행 × 17컬럼
- 한 행 = 화합물 하나와 단백질 하나의 조합(compound-protein pair)
- 컬럼 구성: 화합물 특성(`molecular_weight`, `logp`, `h_bond_donors`, `h_bond_acceptors`, `rotatable_bonds`, `polar_surface_area`, `compound_clogp`), 단백질 특성(`protein_length`, `protein_pi`, `hydrophobicity`), 결합 관련(`binding_site_size`, `mw_ratio`, `logp_pi_interaction`), 식별자(`compound_id`, `protein_id`), 타깃(`binding_affinity`), 파생 이진 컬럼(`active`)
- `compound_id`: 2,000개 전부 고유. `protein_id`: 400개 그룹, 그룹당 평균 5건(표준편차 2.1, 최소 1 / 최대 11)

## 3. 전처리

| 처리 | 내용 | 근거 |
|---|---|---|
| 컬럼 제외 | `active` | `binding_affinity`를 임곗값 7.0으로 이진화한 완전 파생 컬럼(active=0 최댓값 6.996 / active=1 최솟값 7.002) — target leakage |
| 컬럼 제외 | `compound_id` | 완전 고유 식별자, 피처로서 정보 없음 |
| 결측 행 삭제 | `logp`/`polar_surface_area`/`hydrophobicity` 관련 174행 | 각 컬럼 3.0%씩 결측이나 서로 거의 겹치지 않아 합집합 기준 174행(8.7%). 표본이 충분해 삭제 |
| 이상값 재분류 후 삭제 | `polar_surface_area` = -24.65 (`CID_01548`) 1건 | 표면적은 정의상 0 이상 — 물리적으로 불가능한 값. 결측 처리 후 삭제 |
| 중복 확인 | 완전 중복행 0건, `compound_id`+`protein_id` 조합 중복 0건 | 조치 불필요 |
| 데이터 타입 확인 | 위장 결측 문자열 0건, 쉼표/단위 섞인 숫자 0건, ID 형식 100% 일관 | 조치 불필요 |

**행 수**: 2,000 → 1,826(결측 삭제) → **1,825**(이상값 삭제), 15컬럼. `data/processed/drug_discovery_virtual_screening_processed.csv`.

## 4. EDA

**`binding_affinity` 분포** (n=1,825): 평균 6.522, 중앙값 6.476, 표준편차 1.205, 최솟값 1.990, 최댓값 15.040, 오른쪽으로 치우침(skew 0.678). IQR 기준 이상값 58건(3.2%), 상·하단에 대칭적으로 분포.

**상관관계** (`binding_affinity` 기준, 절대값 순)

| feature | r |
|---|---|
| `logp_pi_interaction` | +0.751 |
| `logp` | +0.607 |
| `protein_pi` | +0.295 |
| `mw_ratio` | -0.071 |
| 나머지 9개 | \|r\| < 0.06 |

**다중공선성**: `logp` ↔ `logp_pi_interaction`(r=0.81), `protein_length` ↔ `mw_ratio`(r=-0.82) — 각 쌍은 한쪽이 다른 쪽의 파생값이다.

**분포 특이사항**: `mw_ratio`가 가장 비정상(skew 1.62, 이상값 6.79%로 최다). `protein_length`는 종모양이 아닌 균일분포. `logp` vs `binding_affinity` 산점도에서 주 군집 아래로 별도 띠가 관찰돼 `protein_pi`와의 상호작용 가능성을 초기에 의심했다(7장에서 검증).

## 5. Baseline 모델링

**피처 선정 기준** (해석 가능성 우선): ID 컬럼 제외, target leakage 제외, 결측 과다 변수 제외, 상관 0.8 이상 다중공선성 쌍은 원본만 채택.

- **채택 7개**: `logp`, `protein_pi`, `molecular_weight`, `protein_length`, `hydrophobicity`, `rotatable_bonds`, `polar_surface_area`
- **제외**: `compound_id`/`protein_id`(ID·고카디널리티), `active`(누수), `logp_pi_interaction`/`mw_ratio`(다중공선성), `h_bond_donors`/`h_bond_acceptors`/`binding_site_size`/`compound_clogp`(타깃 상관 \|r\|<0.01)

**학습**: 무작위 8:2 분할(`random_state=42`, train 1,460 / test 365) → `LinearRegression`

| feature | coefficient |
|---|---|
| `logp` | +0.4411 |
| `hydrophobicity` | +0.4184 |
| `protein_pi` | +0.2789 |
| `rotatable_bonds` | +0.0089 |
| `polar_surface_area` | +0.0008 |
| `protein_length` | +0.0002 |
| `molecular_weight` | -0.0001 |
| intercept | 2.6515 |

피처 스케일이 서로 달라(`protein_length` 200~1500 vs `hydrophobicity` 0.3~1.0) 계수 크기를 그대로 중요도로 비교할 수 없다(미표준화).

**성능**: test R² 0.4275, MAE 0.4714, RMSE 0.9916 (train R² 0.4666, MAE 0.4337, RMSE 0.8591). 항상 평균(6.52)만 예측하는 더미 대비 test MAE 약 47%, RMSE 약 24% 낮음.

## 6. 모델 비교

**트리 모델 포함 1차 비교**

| 모델 | train R² | test R² | test MAE | test RMSE |
|---|---|---|---|---|
| Linear Regression (baseline) | 0.4666 | 0.4275 | 0.4714 | 0.9916 |
| RandomForest(`max_depth=6, min_samples_leaf=8`) | 0.6214 | 0.4768 | 0.4603 | 0.9479 |
| XGBoost(`max_depth=3, learning_rate=0.05`) | 0.7585 | 0.4717 | 0.4724 | 0.9525 |

**5-fold 교차검증** (1,825행 전체, `KFold(shuffle=True, random_state=42)`)

| 모델 | R² (mean ± std) |
|---|---|
| Linear Regression | 0.4554 ± 0.1036 |
| RandomForest | 0.4887 ± 0.0816 |
| XGBoost | 0.4915 ± 0.0900 |

판정 규칙(평균 차이 vs 표준편차)상 **세 모델 모두 통계적으로 유의미한 차이 없음** — 단일 분할에서의 트리 모델 우위는 fold 간 변동성(std 0.08~0.10) 범위 안에 있었다. 이에 따라 피처 구성과 정규화 방식을 바꿔 개선 여지를 재검토했다.

**추가 실험 3건** (동일 표본 재현, 이후 교차검증으로 재검증)

| 실험 | test R² | 5-fold CV R² (mean±std) |
|---|---|---|
| `logp`×`protein_pi` 상호작용항 추가 | 0.4469 | 0.4647 ± 0.0454 |
| `protein_id` 그룹 분할(`GroupShuffleSplit`) | 0.4470 | — |
| Ridge 7피처(StandardScaler) | 0.4274 | 0.4552 ± 0.0577 |
| **Ridge 9피처(+`logp_pi_interaction`, +`mw_ratio`)** | **0.5533** | **0.5714 ± 0.0681** |

- **상호작용항**: 교차검증 기준 baseline 대비 차이(0.0095)가 표준편차(0.0577)보다 작아 유의미하지 않음. 7장 residual analysis의 상호작용 기각 결론과 일치.
- **그룹 분할**: `protein_id`는 400개 단백질에 화합물이 평균 5개씩 묶인 구조다. 무작위 분할은 같은 단백질이 train/test에 동시에 들어갈 수 있어 "이미 학습에 등장한 단백질에 새 화합물이 왔을 때"에 가까운 평가이고, `protein_id` 기준 그룹 분할은 train/test의 단백질을 완전히 분리(겹치는 protein_id 0개 확인)해 "한 번도 보지 못한 새 단백질"에 대한 일반화를 측정한다. 두 방식의 성능(0.447 vs 0.427~0.455)이 비슷해 — **이 모델은 새 화합물이든 새 단백질이든 비슷한 수준으로 일반화된다.**
- **Ridge 7 vs 9피처**: 7피처는 규제를 걸어도 baseline과 성능이 거의 동일했다. 반면 다중공선성 때문에 제외했던 `logp_pi_interaction`/`mw_ratio`를 다시 포함한 9피처 Ridge는 R²가 크게 상승했고(0.4552→0.5714), baseline과의 차이(0.1162)가 표준편차(0.0681)보다 커 통계적으로 유의미했다.

**해석용 모델과 예측용 모델**: 5장에서 다중공선성 쌍의 원본만 남긴 것은 계수 하나하나를 안정적으로 해석하기 위한 선택이었다(해석용 모델). Ridge는 계수 크기에 벌점을 부여해 상관 높은 피처를 함께 넣어도 계수가 불안정해지지 않는다는 특성이 있어, 이를 활용해 제외했던 두 피처를 되살린 결과가 9피처 Ridge다(예측용 모델). 두 결과는 상충하지 않는다 — "무엇이 왜 영향을 주는가"를 볼 때는 7피처 모델을, "얼마나 정확히 맞히는가"가 목적일 때는 9피처 Ridge를 쓴다. 이번 목적(스크리닝 성능)에 맞춰 **9피처 Ridge를 최종 모델로 채택**했다. 표준화 계수 기준 `logp_pi_interaction`(+0.933)이 압도적 1위이며, 원본 `logp`(-0.024)/`protein_pi`(-0.022)의 기여는 거의 사라진다 — 데이터셋에 이미 있던 `logp_pi_interaction`은 직접 계산한 `logp`×`protein_pi`(상호작용항 실험, 개선 미미)보다 강한 신호를 담고 있었다는 뜻이다.

**최종 모델**: Ridge(9피처, StandardScaler, alpha=1.0)를 train+test 전체(1,825행)로 재학습해 `output/day6/models/final_model.pkl`로 저장했다. 일반화 성능은 재학습 자체의 적합도가 아니라 앞선 교차검증 결과(**R² 0.5714 ± 0.0681**)를 기준으로 삼는다.

## 7. Residual Analysis

**핵심 메시지: 모델은 결합력이 특히 높은 극단값을 체계적으로 과소예측한다.** 이 편향은 baseline과 최종 모델 모두에서 확인됐고, 이상값은 삭제하지 않고 진단만 했다.

| | Baseline (7피처 LR) | 최종 모델 (Ridge 9피처)† |
|---|---|---|
| corr(residual, 실제값) | **0.773** | **0.736** |
| 최대 오차 사례(실제 15.04) | residual +6.96 | residual +6.46 (약 7% 개선) |
| \|residual\| 최고 상관 feature | `hydrophobicity`(r=0.179) | `hydrophobicity`(r=0.159) |

† held-out 평가를 위해 train만으로 학습한 동일 구조 모델(`ridge9.pkl`)을 기존 test set(365행)에 적용.

- **Baseline**: 오차 최대 10건 전부 타깃 극단값에서 발생 — 실제값 10.4~15.0(4건)은 전부 과소예측, 2.4~2.9(5건)은 전부 과대예측. \|residual\|과 `hydrophobicity`의 상위 4분위(0.715~0.869) 평균 \|residual\|은 0.80으로 나머지 구간(0.30~0.40)의 약 2배. \|residual\|과 `logp×protein_pi`의 상관은 0.034로, 4장에서 의심한 상호작용 구조는 오차를 설명하지 않음(기각).
- **최종 모델**: 실제값 구간별 평균 residual은 하위 10% -0.815(과대예측) / 중간 80% +0.013(거의 무편향) / **상위 10% +0.957(여전히 뚜렷한 과소예측)**.

corr(residual, 실제값)이 0.773→0.736으로만 낮아진 것에서 보듯, **R²가 크게 개선된 뒤에도 이 편향은 거의 그대로 남았다.** 모델 구조를 바꾸는 것만으로는 해결되지 않을 가능성이 높으며, 현재 피처가 담고 있는 정보 자체의 한계로 보인다.

## 8. 한계

- 최종 모델도 `binding_affinity` 변동의 약 43%(1-0.5714)를 설명하지 못한다.
- 스크리닝에서 가장 중요한 고결합력 구간에서 가장 부정확하다(7장) — 유망 후보를 놓칠 위험이 있다.
- 교차검증에서 fold 간 R² 변동폭이 크다(baseline 기준 0.35~0.65) — 단일 분할 비교만으로는 모델 순위를 단정하기 어렵다.
- Ridge의 `alpha`는 튜닝하지 않고 기본값(1.0)을 사용했다.
- `hydrophobicity` 상위 구간에서 오차가 커지는 패턴의 원인은 특정하지 못했다.
- 검증은 모두 이 데이터셋 내부 분할(무작위/그룹) 기준이며, 외부 데이터로는 검증하지 않았다.

## 9. 최종 결론 / Next Steps

**최종 모델로 Ridge 9피처(StandardScaler, alpha=1.0)를 선택한 이유**:
1. 7피처 baseline 대비 교차검증 R²가 0.4552→0.5714로 유의미하게 높다(차이 0.1162 > std 0.0681).
2. RandomForest·XGBoost와도 동등하거나 나은 수준이며, fold 간 변동성(std 0.068)이 더 작아 안정적이다.
3. `protein_id` 그룹 분할로도 성능이 유지돼, 새로운 단백질에 대한 일반화도 근거 있게 확인됐다.
4. 다만 고결합력 극단값 과소예측 편향은 해소되지 않아, 스크리닝 실사용에는 이 한계를 인지한 채로 참고 자료로만 써야 한다.

**남은 과제**:
1. Ridge `alpha` 하이퍼파라미터 튜닝(`GridSearchCV`)
2. 극단값 편향 자체를 겨냥한 접근(분위수 회귀, 가중 손실 등) 검토
3. 외부 데이터로 최종 모델 검증
