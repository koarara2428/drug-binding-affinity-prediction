# 작업 흐름 — drug_discovery_virtual_screening

대상 데이터: `data/raw/drug_discovery_virtual_screening.csv`
목표: `binding_affinity`(결합 친화도)를 `molecular_weight`, `logp`, `protein_pi` 등으로 예측하는 회귀분석. 성공 기준은 "실험 전 컴퓨터 스크리닝 목적의 초보 연구원에게 쓸모 있는 수준".
진행 절차: `docs/reference/SOP.md`의 STAGE ①~⑤ 순서를 따른다.

## 현재 상태

- **완료**: STAGE ① 원본 탐색, STAGE ② 전처리(컬럼 제외 + 결측 삭제 + PSA 이상값 삭제) **최종 완료** — `data/processed/drug_discovery_virtual_screening_processed.csv`(1,825행 × 15컬럼). STAGE ③ 시각화(상관 히트맵/분포/산점도, 전체 feature 히스토그램·boxplot, 최종본 기준 상관분석) 완료. Train/test 분할 방식은 그룹(`protein_id`) 분할 대신 SOP 기본값(무작위 8:2, random_state=42)으로 진행하기로 결정
- **다음 할 일**: 사용자가 STAGE ④(회귀분석) 착수는 별도 지시가 있을 때까지 대기 요청 — 지시 있을 때 진행
- **보류**: 없음 (`polar_surface_area` 이상값 처리 완료로 해소됨)
- **막힌 것**: 없음
- **결론(잠정)**: `binding_affinity` 회귀는 가능. 단 `active`는 타깃의 이진화 파생값이라 반드시 제외, `mw_ratio`/`logp_pi_interaction`은 다른 컬럼의 파생값이라 다중공선성 주의 필요. `logp` 산점도에서 이중 띠 패턴 발견 → `logp`×`protein_pi` 상호작용 가능성, STAGE ④에서 트리 계열 비교 필수

## 진행 내역

### 1. 회귀분석 가능 여부 판정
행 수(2,000), 컬럼 수(17), 연속형 타깃 후보, 결측 비율(최대 3%), 타깃 누수 여부 6개 항목 점검.
→ **결론**: `binding_affinity`를 타깃으로 회귀 가능. `active`는 `binding_affinity`를 임곗값 7.0으로 이진화한 완전 누수 컬럼(active=0 최댓값 6.996 / active=1 최솟값 7.002)이라 피처에서 반드시 제외.
사용 키워드: `pandas.read_csv`, `nunique`, `isna().mean()`, `corr`

### 2. 폴더 정리
`output/day6/` 생성 (요청에 따라 `docs/`는 기존 것 유지, `output/`은 `outputs`가 아니라 실제로는 `output` 단수).

### 3. STAGE ① 원본 탐색 — `docs/reference/SOP.md` 1-1~1-3 그대로 적용
- 1-1 기본 스캔: dtype 17개(수치 15 + str 2), 상수 컬럼 없음, 위장 결측 없음
- 1-2 식별자 확인: `compound_id` 완전 유일(중복 0), `protein_id`는 400개 그룹(그룹당 평균 5건) → 식별자 아닌 그룹 변수
- 1-3 사전 진단: 다중공선성 0.8 이상 쌍 2개 — `protein_length`↔`mw_ratio`(0.814), `logp`↔`logp_pi_interaction`(0.807). 표본 2,000행으로 교차검증 필수 기준(test<100행)에는 해당 없음
- SOP에서 이 데이터에 맞지 않는 항목 정리: 위장 결측 탐지, surrogate ID 생성, 소수 카테고리 기준 `stratify` — 모두 해당 없음. 반대로 `active` 누수 문제는 SOP 체크리스트에 없는 새 발견이라 별도 기록.
사용 키워드: `duplicated`, `is_unique`, `value_counts`, `corr().abs()`

### 4. 분석 주제 5개 브레인스토밍
`binding_affinity`, `molecular_weight`, `polar_surface_area`, `protein_pi`, `binding_site_size`를 각각 타깃으로 한 주제 제시. 오늘 안에 끝낼 후보로 `molecular_weight`(쉬움)와 `binding_affinity`(핵심 주제, 보통)를 추천.

### 5. 목표 확정 + 초기 탐색 3가지
목표: `binding_affinity` ← `molecular_weight`, `logp`, `protein_pi`. 대상: 스크리닝 연구원.
- 최고 `binding_affinity`: `CID_00869`↔`PID_421` (15.04)
- `molecular_weight` 영향: 상관계수 -0.0109, 5분위 평균도 6.50~6.55로 평평함 → **영향 거의 없음**
- `protein_pi` 상관관계: 상관계수 0.2956, 5분위 평균이 6.06→7.05로 꾸준히 증가 → **약한 양의 상관**
사용 키워드: `sort_values`, `qcut`, `groupby`, `corr`

### 6. STAGE ③ 시각화 — `docs/reference/SOP.md` 3-1~3-4 적용 (원본 데이터 기준, 수정 없음)
`scripts/01_eda_structure_overview.py` 작성·실행. 결과는 `output/day6/figures/`에 저장(`correlation_heatmap.png`, `histogram_binding_affinity.png`, `scatter_target_vs_features.png`).
- 이상치(IQR, 전체 기준): `mw_ratio` 144건(7.2%)이 가장 많음, `binding_affinity` 66건(3.3%), `rotatable_bonds`/`molecular_weight`/`logp` 각 30~46건. `h_bond_donors`/`acceptors`/`rotatable_bonds`는 값 종류가 적은 정수 카운트라 IQR 판정을 그대로 삭제 근거로 쓰기엔 이름
- 다중공선성 재확인: `logp`↔`logp_pi_interaction`(0.81), `protein_length`↔`mw_ratio`(-0.81) — STAGE ①과 동일
- `binding_affinity` 분포: 6 부근 주 봉우리 + 오른쪽 긴 꼬리(최대 15.04), `active` 임곗값(7.0)이 오른쪽 어깨를 가로지름
- 산점도(선형성): `molecular_weight`는 추세 없음(r=-0.01), `protein_pi`는 약한 양의 상관(r=0.30), **`logp`는 뚜렷한 양의 상관(r=0.60)이지만 주 군집 아래로 별도 대각선 띠가 존재** → `logp`×`protein_pi` 상호작용 가능성(`logp_pi_interaction` 상관 0.75로 `logp` 단독보다 높은 것과 일치)
- Simpson's Paradox 점검: 소수 카테고리 그룹 변수가 없어 해당 없음(`protein_id`는 400개로 과함, `active`는 타깃 파생이라 사용 불가)
사용 키워드: `matplotlib`, `seaborn.heatmap`, `sns.histplot`, `qcut` 없이 `corr`/`quantile` 기반 IQR

### 7. STAGE ② 전처리 실행 — `scripts/02_preprocessing.py`
- 컬럼 삭제: `active`(타깃 누수), `compound_id`(순수 식별자) 제외 → 15컬럼
- 결측치 처리: `logp`/`polar_surface_area`/`hydrophobicity` 결측 행 삭제로 확정(추천안 채택). **컬럼별로는 각 3%였지만 세 컬럼의 결측이 서로 거의 겹치지 않아 실제 삭제된 행은 174행(8.7%)** — 컬럼별 비율보다 실손실이 큼을 확인·기록
- 결과: 2,000행 → 1,826행. `data/processed/drug_discovery_virtual_screening_processed.csv` 저장 후 shape/결측/제외컬럼 `assert` 재검증 통과
사용 키워드: `drop(columns=...)`, `dropna(subset=...)`, `to_csv`, `assert`

### 8. 전처리 전 문제 7항목 재점검 (기존 항목 3·7과 신규 항목 통합, 중복 서술 없이 정리용으로만 추가)
- 결측치: 항목 7 참고(174행/8.7%, 이미 처리됨)
- **중복 데이터(신규 확인)**: 전체 컬럼 기준 완전 중복행 0건, `compound_id`+`protein_id` 조합 중복도 0건 → 조치 불필요
- 이상값: 항목 6 참고(다수 컬럼 IQR 이상값, `binding_affinity`는 스크리닝 목적상 유지 권장). **신규**: `polar_surface_area`에 물리적으로 불가능한 음수값 1건(-24.65, `CID_01548`) 발견 → 결측 재분류 후 삭제 제안, 실행은 확인 대기
- **잘못된 데이터 타입(신규 확인)**: 쉼표/단위 섞인 숫자, object로 잘못 읽힌 수치 컬럼 없음. ID 형식(CID_/PID_) 100% 일관 → 문제 없음
- ID 컬럼: `compound_id`(완전 유일, 이미 제외) / `protein_id`(400개 그룹 변수, 데이터엔 유지하되 모델 피처로는 미사용 제안)
- 미사용 컬럼: `compound_id`(제외 완료), `protein_id`(피처 제외 제안), `mw_ratio`/`logp_pi_interaction`(원본 피처와 다중공선성 0.81 — 모델링 단계에서 원본/파생 중 택1 제안)
- target leakage: `active`만 완전 누수(이미 제외). `mw_ratio`/`logp_pi_interaction`은 피처 간 파생이라 누수 아님(다중공선성과는 별개 이슈)
- **행 개수**: 원본 2,000행 → 결측 삭제 후 1,826행 → PSA 이상값까지 삭제 후 1,825행 (진행 내역 9번에서 최종 실행)

### 9. 전체 feature 분포 확인 + PSA 이상값 최종 처리 + 상관분석 — `scripts/04_eda_all_numeric_features.py`, `02_preprocessing.py`(갱신), `05_correlation_analysis.py`
- `binding_affinity` 제외 13개 수치형 feature 히스토그램(그리드) + 이상값 상위 4종(`mw_ratio`/`rotatable_bonds`/`logp`/`molecular_weight`) boxplot 작성. `output/day6/figures/numeric_features_histograms.png`, `numeric_features_boxplots_notable.png`
- 특이 분포: `mw_ratio`(skew 1.62, 이상값 6.79%로 최다, 파생값이라 원본 오류 아님), `protein_length`(종모양이 아닌 균일분포)
- 사용자 지시로 `polar_surface_area` 음수값(-24.65, `CID_01548`) **삭제 확정 실행** — `02_preprocessing.py`에 재분류·삭제 단계 추가 후 재실행. 1,826 → **1,825행**으로 최종 확정
- STAGE ④ 착수는 사용자 지시가 있을 때까지 대기하기로 함 — 대신 최종본(1,825행) 기준으로 상관분석 요청받아 진행:
  - correlation matrix 저장: `output/day6/correlation_matrix.csv`
  - heatmap: `output/day6/figures/correlation_heatmap_final.png`
  - `binding_affinity` 상관 랭킹(절대값 기준): `logp_pi_interaction`(+0.75) > `logp`(+0.61) > `protein_pi`(+0.30) > `mw_ratio`(-0.07) > 나머지 |r|<0.06(사실상 무관)
사용 키워드: `corr()`, `sns.heatmap`, `sort_values(key=np.abs)`, `dropna`, `to_csv`

### 10. Linear Regression 피처 7개 확정 + train/test 분할 — `scripts/06_train_test_split.py`
- 상관분석(9번) 기반으로 피처 7개 제안 후 사용자 확정: `logp`, `protein_pi`, `molecular_weight`, `protein_length`, `hydrophobicity`, `rotatable_bonds`, `polar_surface_area` (+ target `binding_affinity`)
- 다중공선성 쌍 처리: `logp_pi_interaction`(↔`logp` 0.81) 제외, `mw_ratio`(↔`protein_length` -0.82) 제외 — 각 쌍의 원본만 채택
- `data/processed/split/` 하위 폴더 생성 후 무작위 8:2 분할(`random_state=42`, 그룹 분할 없음) 실행: train 1,460행 / test 365행. `train.csv`/`test.csv` 저장, target 평균(6.52 vs 6.52) 거의 동일함을 확인
사용 키워드: `train_test_split`, `to_csv`

### 11. Linear Regression 베이스라인 학습·평가 — `scripts/07_linear_regression.py`
- `train.csv`로 학습 → `test.csv`로 예측. 결과: R²(train 0.467 / test 0.428), MAE(0.434 / 0.471), RMSE(0.859 / 0.992) — train-test 격차 0.04 정도로 과적합 크지 않음
- 회귀계수(절대값 순): `logp`(+0.441) ≈ `hydrophobicity`(+0.418) > `protein_pi`(+0.279) > 나머지 4개는 0에 가까움(`rotatable_bonds` +0.009, `polar_surface_area` +0.0008, `protein_length` +0.0002, `molecular_weight` -0.0001). **주의**: 피처 스케일이 제각각(예: `protein_length` 200~1500 vs `hydrophobicity` 0.3~1.0)이라 계수 크기를 그대로 중요도로 비교하면 안 됨 — `hydrophobicity`가 `logp`와 비슷한 계수를 갖는 건 스케일이 작기 때문일 수 있음(→ Ridge 단계에서 `StandardScaler`로 재확인 필요)
- 실제값 vs 예측값 산점도: 중앙 군집은 y=x에 가깝게 따라가지만 상단 이상값(결합력 높은 후보들)에서 과소예측 경향 뚜렷 — 스크리닝 목적상 가장 중요한 구간에서 정확도가 떨어짐을 확인
- 결과 저장: `output/day6/models/linear_regression_metrics.csv`, `linear_regression_coefficients.csv`, `linear_regression.pkl`(model + feature_columns), `output/day6/figures/linear_regression_actual_vs_pred.png`
- **버그 수정**: matplotlib 한글 라벨이 기본 폰트(DejaVu Sans)에서 깨짐 → `font.family="Malgun Gothic"` 지정으로 해결. 이후 시각화 스크립트에도 적용 필요
사용 키워드: `LinearRegression`, `r2_score`, `mean_absolute_error`, `mean_squared_error`, `joblib.dump`

### 12. Linear Regression residual analysis — `scripts/08_residual_analysis.py`
사용자 요청으로 Ridge보다 먼저 진행(비선형 문제인지 다중공선성 문제인지 구분하기 위해). 이상치는 삭제하지 않고 그대로 진단만 함.
- residual 분포: 평균 -0.0005(거의 무편향)이지만 skew 1.67로 오른쪽 꼬리 김. 대부분 -1~+1 사이에 몰려있지만 일부 +7까지 벌어짐
- 예측값 vs residual: 중앙 군집(예측 5~8)은 -1~+1 사이로 고르게 분포(뚜렷한 이분산성 없음), 대신 몇 개 극단치가 예측값 전 구간에 흩어져 있음
- **오차 최대 10개는 전부 타깃 극단값**: 실제값이 아주 높은 경우(10.4~15.0, 4건) 전부 과소예측, 실제값이 아주 낮은 경우(2.4~2.9, 5건) 전부 과대예측 — 중간 구간에는 큰 오차 없음
- **corr(residual, 실제값) = 0.773**로 매우 강함 → 전형적인 "평균으로 회귀" 편향: 극단적으로 높은/낮은 실제값을 모델이 평균 쪽으로 뭉뚱그려 예측
- |residual|과 상관 가장 높은 feature: `hydrophobicity`(0.179) — 4분위 중 최상위 구간(0.715~0.869)만 평균 |residual| 0.80으로 나머지(0.30~0.40)의 약 2배. `protein_pi`도 완만한 증가 경향(0.44)
- **상호작용 가설 기각**: |residual|과 `logp×protein_pi`(제외했던 상호작용항 proxy) 상관 0.034로 사실상 무관 — "상호작용을 못 잡아서 오차가 크다"는 이전 가설은 이 분석에서 뒷받침되지 않음. 진짜 원인은 상호작용이 아니라 **타깃 극단값에서의 구조적 회귀-평균 편향** + `hydrophobicity` 상위 구간
- 결과 저장: `output/day6/figures/residual_distribution.png`, `predicted_vs_residual.png`, `output/day6/models/linear_regression_top10_errors.csv`
사용 키워드: `joblib.load`, `qcut`, `corr`, residual = 실제값 - 예측값

### 13. RandomForest + XGBoost 학습·평가 + 3모델 비교 — `scripts/09_random_forest.py`, `10_xgboost.py`, `11_model_comparison.py`
사용자 요청으로 Ridge 건너뛰고 바로 진행. 피처·target은 Linear Regression과 동일 7개.
- RandomForest(`n_estimators=300, max_depth=6, min_samples_leaf=8`): R² train 0.621 / test 0.477, OOB 0.487(test와 근접해 신뢰할 만함). train-test 격차 0.145 — SOP 기준(0.1) 초과, 경미한 과적합. leaf당 평균 샘플 34.7건으로 개별 암기는 아님
- XGBoost(`max_depth=3, learning_rate=0.05, reg_alpha=0.1, reg_lambda=1.0, subsample/colsample=0.8`): R² train 0.759 / test 0.472, 격차 0.287로 셋 중 과적합이 가장 큼
- **3모델 test 비교**: RandomForest(R² 0.477, MAE 0.460, RMSE 0.948) ≈ XGBoost(0.472, 0.472, 0.953) > Linear Regression(0.428, 0.471, 0.992). RF/XGB 차이는 미미(사실상 동률), 둘 다 선형회귀보다 소폭 우수
- **핵심 한계 재확인**: RF/XGBoost의 실제값 vs 예측값 산점도에서도 타깃 상단 극단값(10~15)이 여전히 대각선 아래로 처짐 — residual analysis에서 발견한 "평균으로의 회귀" 편향이 트리 계열에서도 완전히 해소되지 않음
- feature importance: 두 모델 모두 `logp`가 압도적 1위(RF 0.72 / XGB 0.41), `protein_pi`가 2위 — Linear Regression 계수 순위와 대체로 일치
- 결과 저장: `output/day6/models/{random_forest,xgboost}_metrics.csv`, `*_feature_importance.csv`, `*.pkl`, `model_comparison.csv`, `output/day6/figures/{random_forest,xgboost}_actual_vs_pred.png`, `model_comparison.png`
사용 키워드: `RandomForestRegressor(oob_score=True)`, `XGBRegressor`, `feature_importances_`

### 14. 5-fold 교차검증 — `scripts/12_cross_validation.py`
train+test(1,825행) 재결합 후 `KFold(5, shuffle=True, random_state=42)`로 재검증. 목적: 단일 분할에서 본 RF(0.477)/XGB(0.472) > Linear(0.428) 순위가 진짜인지 확인.
- 결과: Linear Regression R² 0.455±0.104 / RandomForest 0.489±0.082 / XGBoost 0.492±0.090
- **SOP 4-7 판정 규칙(mean 차이 vs std) 적용 결과, 세 모델 모두 "사실상 동률"** — Linear vs RF 차이 0.033 < std 0.104, RF vs XGB 차이 0.003 < std 0.090. **이전에 단일 분할로 내린 "RandomForest가 가장 낫다"는 결론은 fold 간 변동성(R² std 0.08~0.10, boxplot상 최저 0.35~최고 0.65까지 폭넓게 흩어짐)에 비해 근거가 약했다는 게 드러남 — 정정 필요**
- 결과 저장: `output/day6/models/cross_validation_results.csv`, `output/day6/figures/cross_validation_r2_boxplot.png`
사용 키워드: `KFold`, `cross_validate`, `scoring=["r2","neg_mean_absolute_error","neg_root_mean_squared_error"]`

### 15. STAGE ⑤ 문서화 — 루트 `README.md` 작성
사용자 지정 11개 섹션 구조(Objective~Next steps)로 오늘 수행한 분석만 정리. 실제 계산된 수치만 사용, 과장 없이 작성. 성능 비교 결과(교차검증에서 3모델 동률)와 한계(회귀-평균 편향, `protein_id` 그룹 미반영, Ridge 미실행 등)를 숨기지 않고 포함.

### 16. README "Next steps" 4가지 실행 — `scripts/13_next_steps.py`, `14_next_steps_cv_check.py`, `15_final_model.py`
- 1) 상호작용항(`logp`×`protein_pi`) 추가 LR: R² 0.455→0.465(CV), 차이가 std보다 작아 유의미하지 않음 — 8번 residual analysis의 상호작용 기각 결론과 일치
- 2) `protein_id` GroupShuffleSplit(겹치는 protein_id 0개 확인): R² 0.447, 기존 무작위 분할과 비슷하거나 소폭 나음 — 무작위 분할의 단백질 정보 누수 우려는 확인되지 않음
- 3) Ridge+StandardScaler: 7피처는 일반 LR과 거의 동일(R² 0.455). **9피처(+`logp_pi_interaction`,+`mw_ratio`)는 R² 0.455→0.571로 크게 상승, CV로 유의미함 확인(차이 0.116 > std 0.068)**. 표준화 계수에서 `logp_pi_interaction`(+0.933)이 압도적, 원본 `logp`/`protein_pi` 기여는 거의 사라짐
- 4) **최종 모델 확정: Ridge(9피처, StandardScaler, alpha=1.0)**. train+test 전체(1,825행)로 재학습 후 `output/day6/models/final_model.pkl` 저장. 일반화 성능은 CV 기준 R² 0.571±0.068 사용(재학습 자체 적합도는 참고용일 뿐 과적합 여부 판단 불가라고 스크립트 출력에 명시)
- `README.md`에 "12. Next steps 수행 결과" 섹션 추가, 5장(Feature selection)의 원래 결론이 갱신됐음을 명시
사용 키워드: `GroupShuffleSplit`, `StandardScaler`, `Ridge`, `make_pipeline`, `KFold`

### 17. 최종 모델(Ridge 9피처) residual 재검증 — `scripts/16_final_model_residual_analysis.py`
- **버그 발견·수정**: `ridge9.pkl`(script 13)은 `scaler`와 `model`을 파이프라인이 아니라 따로 저장했는데, 처음에 `model.predict(원본피처)`로 스케일링을 건너뛰어 예측값이 140대로 폭주하는 오류 발생 → `scaler.transform()`을 먼저 적용하도록 수정 후 재실행
- `final_model.pkl`은 train+test 전체로 재학습돼 held-out 데이터가 없으므로, train만으로 학습한 `ridge9.pkl`을 기존 test set(365행)에 적용해 진단
- corr(residual, 실제값) = **0.736** (기존 LR 7피처 0.773 대비 소폭 개선에 그침)
- 실제값 구간별 평균 residual: 하위10% -0.815 / 중간80% +0.013 / **상위10% +0.957** — 극단값 과소예측 편향이 거의 그대로 남아있음 확인
- 오차 최대 사례(실제 15.04): residual 6.46 (기존 LR 6.96 대비 약 7% 개선에 그침)
- \|residual\|과 상관 최고 feature는 이번에도 `hydrophobicity`(0.159) — 8번 항목 결과와 일치
- **결론**: R²는 크게 개선(0.43→0.57)됐지만 스크리닝의 핵심 문제(고결합력 후보 과소예측)는 구조적으로 거의 해결되지 않음. `README.md`에 "13. 최종 모델 residual 재검증" 섹션 추가
- `output/day6/final/`에 `final_model_residual_distribution.png`, `final_model_predicted_vs_residual.png`, `final_model_top10_errors.csv` 추가
사용 키워드: `scaler.transform`, `qcut`, `corr`

### 18. README.md 최종 보고서 형태로 재구성 (내용·수치 변경 없음, 구조만 재정리)
- 최상단에 Executive Summary 추가(최종 모델, CV R² 0.5714±0.0681, baseline 분리 명시)
- 장 순서를 문제정의→데이터셋→전처리→EDA→Baseline 모델링→모델 비교→Residual Analysis→한계→최종결론/Next Steps로 재배열. 기존 11~13장(Next steps 계획/실행/재검증)을 6장(모델 비교)과 7장(Residual Analysis)에 통합
- "해석용 모델 vs 예측용 모델" 프레이밍으로 `logp_pi_interaction`/`mw_ratio`를 뺐다가 Ridge에서 다시 넣은 논리를 명시적으로 설명
- `protein_id` 그룹 구조: 무작위 분할(이미 본 단백질+새 화합물) vs 그룹 분할(완전히 새 단백질)의 의미 차이와, 두 결과가 비슷했다는 한 줄 결론 추가
- 중복 서술 정리, baseline/최종모델 용어 통일
- `output/day6/final/INDEX.md`의 "관련 README" 장 번호를 새 구조에 맞게 갱신
사용 키워드: 없음(문서 편집만)

### 19. README.md를 다른 저장소(bike-sharing-demand-analysis) 형식으로 재구성 + GitHub 업로드
- 사용자의 기존 GitHub 저장소(`koarara2428/bike-sharing-demand-analysis`)를 `gh api`로 조회해 형식을 참고: 상단 소개+목적 문단 → 목차 → 데이터셋(컬럼표) → 프로젝트 구조(트리) → 분석 워크플로우(단계별 설명+**figure 인라인 삽입**) → 모델 성능 비교표 → 핵심 인사이트 → 한계 및 유의사항 → 재현 방법(설치+실행순서+폰트 참고) → Claude Code 활용 방식 → 데이터 출처
- 이 형식에 맞춰 `README.md` 전체 재작성(수치·사실관계는 이전 버전과 동일, Executive Summary는 별도 절 대신 도입부 문단에 통합). Figure 4개 인라인 삽입: `correlation_heatmap_final.png`, `binding_affinity_histogram.png`, `model_comparison.png`, `final_model_predicted_vs_residual.png`
- `requirements.txt` 신규 생성(스크립트 실제 import 기준: pandas/numpy/matplotlib/seaborn/scikit-learn/xgboost/joblib)
- `output/day6/final/INDEX.md`의 "관련 README" 열을 장 번호 대신 새 구조의 섹션명으로 갱신
- **GitHub에 업로드**: `git init` → `.gitignore` 생성(캐시/OS 파일만 제외, 데이터·모델·그림 전부 포함, 총 5.4MB) → 커밋 → `gh repo create koarara2428/drug-binding-affinity-prediction --public --push`로 저장소 생성 및 최초 푸시(84개 파일). 이후 README 재구성은 아직 재푸시 전(사용자 확인 필요)
사용 키워드: `gh api`, `gh repo create`, `git init/add/commit`, `<img>` 대신 마크다운 `![]()` 삽입

### 20. README "데이터 출처" 정정
사용자가 실제 출처를 알려줌: [Kaggle - Drug Discovery Virtual Screening Dataset](https://www.kaggle.com/datasets/shahriarkabir/drug-discovery-virtual-screening-dataset) (Shahriar Kabir 게시). "출처 불명"으로 적었던 문장을 실제 출처 링크로 교체.

## 현재 상태 (최종)

STAGE ①~⑤ + README Next steps 4가지 + 최종 모델 residual 재검증까지 완료. **최종 모델은 Ridge(9피처, StandardScaler) — CV R² 0.571±0.068**, `output/day6/models/final_model.pkl`. **단, 극단값(고결합력 후보) 과소예측 편향은 구조적으로 남아있어 스크리닝 목적에는 참고용 이상으로 쓰기 어렵다는 게 최종 결론.** 산출물: `data/processed/`(전처리본 + train/test), `output/day6/`(그림·모델·지표, `final/`에 핵심만 curation), 루트 `README.md`(Executive Summary + 1~9장, 최종 보고서 형태), `docs/WORKFLOW.md`(과정 로그, 이 파일).
