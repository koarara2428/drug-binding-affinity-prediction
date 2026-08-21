<!-- title: Claude Code로 수행한 신약 탐색 데이터 분석 -->
# Claude Code로 수행한 신약 탐색 데이터 분석 (AI-Assisted Drug Discovery Data Analysis with Claude Code)

> 가상 스크리닝(Virtual Screening) 데이터셋 기반 결합 친화도(Binding Affinity) 예측 및 모델 검증

이 저장소는 데이터 품질 검토, 탐색적 분석(EDA), 모델링, 교차검증, residual 분석까지 이어지는 바이오 데이터 분석 워크플로우 전체를, [Claude Code](https://claude.com/product/claude-code)를 코딩·분석 도구로 활용해 수행한 과정을 보여줍니다. 데이터 품질 판단, 통계적 의사결정, 결과 해석은 모두 연구자가 직접 지시하고 검증했습니다. 이 저장소가 전달하려는 메시지는 "Claude Code가 모델을 만들었다"가 아니라, **"Claude Code를 활용해 분석을 빠르게 수행하면서도, 연구 판단·검증·해석은 연구자가 직접 통제했다"**는 것입니다.

## 목차

- [프로젝트 요약](#프로젝트-요약)
- [수행한 역량](#수행한-역량)
- [핵심 결과](#핵심-결과)
- [분석 워크플로우](#분석-워크플로우)
- [핵심 발견](#핵심-발견)
- [실험 기록](#실험-기록)
- [모델 비교](#모델-비교)
- [핵심 한계](#핵심-한계)
- [Claude Code 활용 방식](#claude-code-활용-방식)
- [연구자 vs AI](#연구자-vs-ai)
- [바이오 연구 맥락](#바이오-연구-맥락)
- [재현 방법](#재현-방법)
- [저장소 구조](#저장소-구조)
- [다음 실험](#다음-실험)
- [상세 문서](#상세-문서)

## 프로젝트 요약

- **Dataset:** 화합물-단백질 조합 2,000건 (품질 검토 후 1,825건)
- **Task:** 결합 친화도(binding affinity) 회귀
- **Workflow:** 데이터 품질 검토 → EDA → 모델링 → 교차검증 → Residual 분석
- **Best model:** Ridge 회귀(9피처, `StandardScaler`)
- **CV R²:** 0.5714 ± 0.0681 (baseline 7피처 Linear Regression: 0.4552 ± 0.0577)
- **Key finding:** 결합력이 높은 화합물을 baseline과 최종 모델 모두에서 체계적으로 과소예측
- **Main purpose:** Claude Code를 활용한, 연구자 검증 기반의 end-to-end 바이오 데이터 분석 워크플로우 시연

## 수행한 역량

- **Claude Code 기반 end-to-end 데이터 분석** — 원본 파일 확인부터 마지막 residual 검증까지 전 단계를 Claude Code를 통해 수행했으며, 한 번에 자동 실행한 것이 아니라 단계마다 방향을 지시했습니다.
- **데이터 품질 검토** — 모델링 전에 원본 파일의 중복, 위장된 결측값, 타입 불일치를 먼저 확인했습니다(`scripts/01_eda_structure_overview.py`).
- **결측치/이상값 처리** — 행의 8.7%에서 결측치를, 물리적으로 불가능한 값(`polar_surface_area` = -24.65) 1건을 발견해 각각에 대해 명시적으로 유지/삭제를 판단했습니다(`scripts/02_preprocessing.py`).
- **데이터 누수 탐지** — `active` 컬럼이 타깃을 그대로 이진화한 값임을 발견해 모델링 전에 제거했습니다 — [핵심 한계](#핵심-한계), [핵심 발견](#핵심-발견) 참고.
- **탐색적 데이터 분석** — 피처 선정 전에 모든 피처의 타깃 상관관계 순위를 매기고 다중공선성 쌍을 표시했습니다(`scripts/03`~`05`).
- **피처 엔지니어링** — `logp`×`protein_pi` 상호작용항을 직접 만들어 검증했고, 다중공선성 때문에 제외했던 피처를 Ridge 정규화로 다시 포함시키는 방법을 평가했습니다(`scripts/13`~`15`).
- **회귀 모델 비교** — Linear Regression, RandomForest, XGBoost를 동일한 분할 조건에서 비교했습니다(`scripts/07`, `09`~`11`).
- **교차검증** — 단일 분할 기준 모델 순위를 5-fold 교차검증으로 재검증했고, 그 결과 최초 결론(어떤 모델이 더 나은가)이 뒤집혔습니다(`scripts/12`, `14`).
- **Residual 분석** — 평균 오차뿐 아니라 *어디서* 틀리는지를 baseline과 최종 모델 모두에서 확인했습니다(`scripts/08`, `16`).
- **모델 한계 평가** — R² 수치만 보고하지 않고, residual 분석 결과를 근거로 최종 모델의 실전 활용 범위를 제한했습니다.
- **재현 가능한 워크플로우 문서화** — 모든 스크립트를 실행 순서대로 번호를 매겼고, 모든 의사결정을 [`docs/WORKFLOW.md`](docs/WORKFLOW.md)에 기록했습니다.

## 핵심 결과

```text
Ridge Regression

7 features
CV R² ≈ 0.455

        ↓ 피처 재검토

9 features
CV R² ≈ 0.571
```

두 피처(`logp_pi_interaction`, `mw_ratio`)는 다중공선성 우려(다른 피처와 상관 0.8 이상) 때문에 baseline 모델에서 처음엔 제외했습니다 — 첫 모델의 계수를 개별적으로 해석 가능하게 유지하기 위한 선택이었습니다. 이후 계수 크기를 규제해 수동 제외 없이도 상관 높은 피처를 함께 쓸 수 있는 Ridge 회귀를 대안으로 검증했습니다. 제외했던 두 피처를 Ridge 모델에 다시 포함시키자 교차검증 R²가 0.4552 ± 0.0577에서 0.5714 ± 0.0681로 상승했고, 이 차이는 두 모델의 fold 간 표준편차보다 컸습니다.

## 분석 워크플로우

```text
원본 데이터 (data/raw/)
   ↓
데이터 품질 검토                → scripts/01_eda_structure_overview.py
   ↓
전처리                          → scripts/02_preprocessing.py
   ↓
EDA                             → scripts/03_eda_binding_affinity.py
                                   scripts/04_eda_all_numeric_features.py
                                   scripts/05_correlation_analysis.py
   ↓
Baseline 모델링                 → scripts/06_train_test_split.py
                                   scripts/07_linear_regression.py
   ↓
모델 비교                       → scripts/09_random_forest.py
                                   scripts/10_xgboost.py
                                   scripts/11_model_comparison.py
   ↓
교차검증                        → scripts/12_cross_validation.py
   ↓
피처/모델 최적화                 → scripts/13_next_steps.py
                                   scripts/14_next_steps_cv_check.py
                                   scripts/15_final_model.py
   ↓
Residual 분석                    → scripts/08_residual_analysis.py (baseline)
                                    scripts/16_final_model_residual_analysis.py (최종 모델)
   ↓
실전 활용성 평가                 → 이 README(핵심 한계), docs/next_experiment.md
```

## 핵심 발견

### 1. 데이터 품질

`logp`, `polar_surface_area`, `hydrophobicity`의 결측치(각 3.0%)는 서로 거의 겹치지 않아 합쳐서 174행(8.7%)에 영향을 줬습니다. 별도 점검에서 물리적으로 불가능한 값도 1건 발견했습니다 — `polar_surface_area` = -24.65(표면적은 음수가 될 수 없음). 둘 다 해당 행을 삭제해 처리했습니다(2,000 → 1,825행). 중복 행이나 위장된 결측 문자열은 발견되지 않았습니다.

### 2. 데이터 누수

`active` 컬럼은 타깃을 그대로 이진화한 값이었습니다: `binding_affinity > 7.0`이면 정확히 `active = 1`(active=0 최댓값 6.996, active=1 최솟값 7.002 — 겹침 없음). 이걸 피처로 쓰면 모델이 답을 미리 보게 되므로 모든 모델링에서 제외했습니다.

### 3. 모델 성능

단일 8:2 분할에서는 RandomForest·XGBoost가 Linear Regression보다 나아 보였습니다. 5-fold 교차검증 결과 이 차이는 신뢰할 수 없었습니다 — fold 간 변동성을 고려하면 세 모델의 R² 분포가 겹쳤습니다. 통계적으로 유의미한 개선을 만든 유일한 변화는 다중공선성 피처 2개를 Ridge 모델에 다시 포함시킨 것이었습니다([핵심 결과](#핵심-결과) 참고).

### 4. 모델의 한계

held-out test set에서 residual 분석 결과 `corr(residual, 실제값)`이 baseline 모델은 0.77, 최종 모델은 0.74였습니다 — 두 경우 모두 결합력이 높은 화합물은 과소예측하고 낮은 화합물은 과대예측합니다. Ridge 모델의 전체 R² 개선은 이 편향을 해소하지 못했고, 이는 실전 활용에 제약이 됩니다 — [핵심 한계](#핵심-한계) 참고.

## 실험 기록

| 실험 | 접근법 | 결과 | 결정 |
|---|---|---|---|
| Baseline | Linear Regression, 7피처 | CV R² = 0.4552 ± 0.0577 | Baseline 확정 |
| 모델 비교 | RandomForest | CV R² = 0.4887 ± 0.0816 | 비교 완료 — baseline 대비 통계적으로 유의미한 개선 없음 |
| 모델 비교 | XGBoost | CV R² = 0.4915 ± 0.0900 | 비교 완료 — baseline·RandomForest 대비 통계적으로 유의미한 개선 없음 |
| 피처 엔지니어링 | `logp`×`protein_pi` 상호작용항 | CV R² = 0.4647 ± 0.0454 | 배운 점: 상호작용항 단독으로는 유의미한 개선 없음 — 최종 모델에서 제외 |
| 데이터 분할 전략 | `protein_id` 기준 그룹 분할 | test R² = 0.4470, 무작위 분할과 유사 | 배운 점: 무작위 분할의 단백질 중복이 성능을 부풀렸다는 근거는 없음 |
| 정규화 | Ridge, 동일 7피처 | CV R² = 0.4552 ± 0.0577 (변화 없음) | 이 피처 조합에서는 기본값(`alpha=1.0`) 규제 효과 없음 |
| 정규화 + 피처 재포함 | Ridge, 9피처(+`logp_pi_interaction`, +`mw_ratio`) | CV R² = 0.5714 ± 0.0681 | 최종 모델로 선정 |
| Residual 분석 | 최종 모델, held-out test set | 고결합력 화합물 여전히 과소예측(`corr(residual, 실제값)` = 0.7358) | 배운 점: R² 개선이 이 편향을 해결하지 못함 — 활용 범위 제한 |

이 표는 **반복적인 연구 워크플로우**를 보여줍니다: 처음부터 정해진 계획이 아니라, 각 실험의 결과가 다음 실험의 방향을 결정했습니다.

## 모델 비교

| 모델 | CV R² (mean ± std) | Test R² (단일 분할) | Test MAE | Test RMSE |
|---|---|---|---|---|
| Linear Regression (baseline, 7피처) | 0.4554 ± 0.1036 | 0.4275 | 0.4714 | 0.9916 |
| RandomForest | 0.4887 ± 0.0816 | 0.4768 | 0.4603 | 0.9479 |
| XGBoost | 0.4915 ± 0.0900 | 0.4717 | 0.4724 | 0.9525 |
| Ridge (7피처) | 0.4552 ± 0.0577 | 0.4274 | 0.4716 | 0.9916 |
| **Ridge (9피처, 최종 모델)** | **0.5714 ± 0.0681** | **0.5533** | **0.3700** | **0.8759** |

*참고: Linear Regression/RandomForest/XGBoost 행은 한 번의 5-fold 교차검증(`scripts/12`)에서, Ridge 실험 행들은 동일한 1,825행 데이터를 다른 행 순서로 재정렬해 실행한 두 번째 5-fold 교차검증(`scripts/14`)에서 나온 값입니다. baseline의 CV 평균이 두 표에서 미세하게 다른 것(0.4554 vs 0.4552)은 이 때문이며, 계산 오류가 아니라 둘 다 같은 모델에 대한 타당한 추정치입니다.*

전체 결과 파일: [`output/day6/final/cross_validation_results.csv`](output/day6/final/cross_validation_results.csv), [`output/day6/final/next_steps_cv_results.csv`](output/day6/final/next_steps_cv_results.csv), [`output/day6/final/model_comparison.png`](output/day6/final/model_comparison.png)

![모델 비교 차트](output/day6/figures/model_comparison.png)

## 핵심 한계

![최종 모델: 예측값 vs Residual, held-out test set](output/day6/figures/final_model_predicted_vs_residual.png)

위 그림(최종 모델, held-out test set)이 핵심 한계를 가장 잘 보여줍니다: 대부분의 점은 residual 0 근처에 촘촘히 몰려 있지만, 뚜렷하게 큰 오차 몇 건이 남아있습니다 — 이는 실제 결합력이 예측보다 훨씬 높거나 낮았던 화합물들이며, [핵심 발견](#핵심-발견)에서 보고한 `corr(residual, 실제값)` = 0.74와 정확히 같은 패턴입니다.

- **데이터셋이 실험 데이터가 아니라 시뮬레이션 데이터입니다.** `data/raw/drug_discovery_virtual_screening.csv`는 가상 스크리닝(virtual-screening) 데이터셋으로([데이터 출처](#상세-문서) 참고), 실험실 측정값이 아니라 컴퓨터로 생성된 값입니다.
- **외부 검증이 없습니다.** 보고된 모든 지표는 이 하나의 데이터셋을 (무작위/그룹 기준으로) 나눈 것이며, 독립적인 다른 데이터셋으로는 평가하지 않았습니다.
- **피처가 기본적인 분자/단백질 서술자로 한정됩니다** (분자량, logP, 단백질 길이, 등전점 등) — 구조 정보나 도킹 기반 피처는 사용하지 않았습니다.
- **결합력이 높은 화합물을 체계적으로 과소예측합니다** (`corr(residual, 실제값)` ≈ 0.74~0.77, 테스트한 모든 모델에서 공통). 이 프로젝트의 모든 모델(최종 모델 포함)에서 동일하게 관찰된 실패 패턴입니다.

**따라서 이 모델을 최종 신약 후보 선별이나 실제 신약개발 의사결정에 바로 사용하는 것은 적절하지 않습니다.** 이 모델은 최종 후보 순위를 매기는 도구가 아니라, 명백히 결합력이 낮을 것으로 확신되는 후보를 걸러내는 **1차 필터링 도구**로 자리매김하는 것이 더 적절합니다. 실제 스크리닝 워크플로우에서 이 프로젝트가 차지하는 위치는 [`docs/biological_context.md`](docs/biological_context.md)를, 고결합력 편향 문제를 해결하려면 무엇이 더 필요한지는 [`docs/next_experiment.md`](docs/next_experiment.md)를 참고하세요.

## Claude Code 활용 방식

```text
연구 질문
      ↓
Prompt / Instruction
      ↓
Claude Code
      ↓
코드 생성 / 분석 실행
      ↓
결과 확인
      ↓
연구자 검증
      ↓
다음 분석 방향 결정
```

이 프로젝트의 모든 단계는 이 루프를 따랐습니다: 구체적인 질문을 던지면 Claude Code가 코드를 생성하고 분석을 실행했고, 연구자는 실제 결과물(수치, 표, 그래프)을 직접 확인한 뒤에야 그 결과를 받아들이고 다음 단계로 갈지, 방향을 바꿀지 판단했습니다. 이 프로젝트의 교차검증 단계가 대표적인 예입니다 — 단일 분할 결과를 그대로 받아들이지 않고 재검증했고, 그 결과 결론이 바뀌었습니다.

각 분석 단계별로 이 루프가 어떻게 적용됐는지는 [`docs/CLAUDE_WORKFLOW.md`](docs/CLAUDE_WORKFLOW.md)에 단계별로 정리했습니다.

## 연구자 vs AI

### Claude Code가 담당한 것

- Python 코드 생성
- 데이터 점검
- 시각화
- 모델 구현
- 반복적인 분석(여러 모델·피처 조합에 걸친 재검증)
- 문서화 지원

### 연구자가 판단한 것

- 문제 정의와 성공 기준
- 데이터 품질 기준(무엇을 결측/무효/중복으로 볼 것인가)
- 데이터 누수 여부 판단
- 전처리 전략(삭제 vs 대체, 어느 기준에서)
- 피처 선정 전략(다중공선성 쌍 중 어느 것을 남길지, 언제 정규화로 둘 다 되살릴지)
- 모델 비교 전략(단일 분할을 신뢰하지 않고 교차검증을 요구)
- 검증 전략(무작위 분할과 그룹 분할을 왜 둘 다 확인했는지)
- 생물학적 해석(residual 편향이 스크리닝 용도에 어떤 의미인지)
- 최종 모델 선정
- 실전 활용 가능성 평가(1차 필터 vs 최종 순위 도구)

## 바이오 연구 맥락

결합 친화도 예측은 신약 개발 초기 단계에서, 컴퓨터로 스크리닝한 화합물 중 어떤 것을 실험으로 넘길지 우선순위를 매기는 데 쓰입니다 — 실험을 대체하는 게 아니라 걸러내는 단계입니다. 이 프로젝트는 순수하게 계산 분석이며, 실제 실험은 수행하지 않았습니다. 이 모델의 한계가 실제 워크플로우의 어느 지점에서 활용 가능한지에 대한 자세한 논의는 [`docs/biological_context.md`](docs/biological_context.md)를 참고하세요.

## 재현 방법

### 환경

```bash
pip install -r requirements.txt
```

### 사용 도구·방법

**AI 활용 개발**
- Claude Code

**프로그래밍**
- Python

**데이터 분석**
- pandas
- NumPy
- matplotlib / seaborn

**머신러닝**
- scikit-learn (Linear Regression, Ridge, RandomForest)
- XGBoost

**검증**
- Train/Test 분할 (8:2, `random_state=42`)
- 5-fold 교차검증
- Residual 분석

### 실행 순서

스크립트는 저장소 루트에서 실행합니다. 스크립트 번호(`01`~`16`)가 곧 실제 의존관계 순서입니다.

```bash
python scripts/01_eda_structure_overview.py
python scripts/02_preprocessing.py
python scripts/03_eda_binding_affinity.py
python scripts/04_eda_all_numeric_features.py
python scripts/05_correlation_analysis.py
python scripts/06_train_test_split.py
python scripts/07_linear_regression.py
python scripts/08_residual_analysis.py
python scripts/09_random_forest.py
python scripts/10_xgboost.py
python scripts/11_model_comparison.py
python scripts/12_cross_validation.py
python scripts/13_next_steps.py
python scripts/14_next_steps_cv_check.py
python scripts/15_final_model.py
python scripts/16_final_model_residual_analysis.py
```

모든 모델이 `random_state=42`로 고정돼 있어, 위 순서대로 실행하면 `output/day6/`의 모든 결과가 동일하게 재현됩니다.

> **폰트 참고**: 일부 그래프의 한글 라벨은 Windows 기본 폰트인 `Malgun Gothic`을 사용합니다. macOS/Linux에서는 각 스크립트의 `plt.rcParams["font.family"]` 값을 시스템에 설치된 한글 폰트(예: `AppleGothic`, `NanumGothic`)로 바꿔야 해당 라벨이 정상적으로 표시됩니다.

### 원본 / 전처리 데이터

- `data/raw/` — 원본 파일, 절대 수정하지 않음.
- `data/processed/` — 정제된 데이터셋(`drug_discovery_virtual_screening_processed.csv`, 1,825 × 15)과 모든 단일 분할 평가에 쓴 8:2 train/test 분할(`split/train.csv`, `split/test.csv`).

### 산출물

- `output/day6/figures/` — 생성된 모든 차트.
- `output/day6/models/` — 프로젝트 진행 중 만든 모든 학습 모델(`.pkl`)과 지표 파일.
- `output/day6/final/` — 위 산출물 중 핵심만 추린 것(최종 모델, 주요 비교표, 주요 그림) + 1페이지 결과 요약 [`FINAL_RESULTS.md`](output/day6/final/FINAL_RESULTS.md) + 파일 안내 [`INDEX.md`](output/day6/final/INDEX.md).

## 저장소 구조

```
.
├── README.md
├── CLAUDE.md
├── requirements.txt
├── data/
│   ├── raw/                          # 원본 파일 (수정 금지)
│   └── processed/                    # 정제된 데이터 + train/test 분할
├── scripts/                          # 분석 스크립트, 01-16 (실행 순서 = 의존관계 순서)
├── output/
│   └── day6/
│       ├── figures/                  # 생성된 모든 차트
│       ├── models/                   # 학습된 모델(.pkl) + 지표
│       └── final/                    # 핵심 결과 curation + FINAL_RESULTS.md + INDEX.md
└── docs/
    ├── WORKFLOW.md                   # 연구 의사결정 로그 (시간순)
    ├── CLAUDE_WORKFLOW.md            # Claude Code 활용 방식 (대표 워크플로우)
    ├── biological_context.md         # 실제 스크리닝 워크플로우에서의 위치
    ├── next_experiment.md            # 제안된(아직 미실행) 후속 실험
    ├── problem_definition.md         # 최초 문제 정의
    ├── reference/SOP.md              # 범용 분석 절차 (참고용)
    └── archive/                      # 이번 분석과 무관한 예전 프로젝트 기록 (보관용)
```

## 다음 실험

[핵심 한계](#핵심-한계)에서 확인한 고결합력 과소예측 편향이 가장 중요한 미해결 과제입니다. 제안된(아직 미실행) 후속 작업으로는 이 편향을 직접 겨냥한 가중/분위수 회귀, XGBoost 하이퍼파라미터 튜닝, 더 풍부한 구조 기반 피처, 외부 검증이 있습니다. 어떤 지표로 실제 개선 여부를 측정할지(Spearman 상관, Top-K recall 포함)는 [`docs/next_experiment.md`](docs/next_experiment.md)에 자세히 정리했습니다.

## 상세 문서

| 문서 | 내용 |
|---|---|
| [`docs/WORKFLOW.md`](docs/WORKFLOW.md) | 전체 시간순 연구 의사결정 로그 |
| [`docs/CLAUDE_WORKFLOW.md`](docs/CLAUDE_WORKFLOW.md) | Claude Code 활용 방식, 단계별 대표 워크플로우 |
| [`docs/biological_context.md`](docs/biological_context.md) | 이 문제가 왜 중요한지, 실제 스크리닝 워크플로우에서 계산 분석이 차지하는 위치 |
| [`docs/next_experiment.md`](docs/next_experiment.md) | 제안된 후속 실험 (아직 미실행) |
| [`docs/problem_definition.md`](docs/problem_definition.md) | 최초 문제 정의 |
| [`output/day6/final/FINAL_RESULTS.md`](output/day6/final/FINAL_RESULTS.md) | 1페이지 결과 요약 |
| [`output/day6/final/INDEX.md`](output/day6/final/INDEX.md) | 핵심 결과 파일 안내 |

**데이터 출처**: [Drug Discovery Virtual Screening Dataset](https://www.kaggle.com/datasets/shahriarkabir/drug-discovery-virtual-screening-dataset) (Kaggle, Shahriar Kabir 게시) — 컴퓨터로 생성된 가상 스크리닝 데이터셋이며, 실험으로 측정한 결합력 값이 아닙니다.
