"""README.md 11장 Next steps 1~3 실행
1) logp x protein_pi 상호작용 항 추가 Linear Regression
2) protein_id 기준 GroupShuffleSplit 재검증
3) Ridge(StandardScaler) — 7피처 버전 / 7+logp_pi_interaction+mw_ratio(9피처) 버전
기존 무작위 분할(train.csv/test.csv)과 동일한 표본을 재현하기 위해 동일 random_state로 전체 처리본을 다시 분할한다.
"""
import pandas as pd
import numpy as np
import joblib
from pathlib import Path
from sklearn.model_selection import train_test_split, GroupShuffleSplit
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error

BASE_DIR = Path(__file__).resolve().parent.parent
PROCESSED = BASE_DIR / "data" / "processed" / "drug_discovery_virtual_screening_processed.csv"
MODEL_DIR = BASE_DIR / "output" / "day6" / "models"

FEATURES = [
    "logp", "protein_pi", "molecular_weight",
    "protein_length", "hydrophobicity", "rotatable_bonds", "polar_surface_area",
]
TARGET = "binding_affinity"

def metrics(y_true, y_pred):
    return {
        "R2": r2_score(y_true, y_pred),
        "MAE": mean_absolute_error(y_true, y_pred),
        "RMSE": np.sqrt(mean_squared_error(y_true, y_pred)),
    }

df = pd.read_csv(PROCESSED)  # 1,825행, protein_id/logp_pi_interaction/mw_ratio 포함

# 기존 8:2 무작위 분할과 동일한 행 구성을 재현 (동일 random_state, 동일 입력 순서)
train_df, test_df = train_test_split(df, test_size=0.2, random_state=42)
print(f"재현된 분할: train {len(train_df)} / test {len(test_df)}")

results = []

# ---------- 1) 상호작용 항 추가 LR ----------
X_train = train_df[FEATURES].copy()
X_test = test_df[FEATURES].copy()
X_train["logp_x_protein_pi"] = train_df["logp"] * train_df["protein_pi"]
X_test["logp_x_protein_pi"] = test_df["logp"] * test_df["protein_pi"]

lr_inter = LinearRegression().fit(X_train, train_df[TARGET])
pred = lr_inter.predict(X_test)
m = metrics(test_df[TARGET], pred)
results.append({"experiment": "1) LR + logp*protein_pi 상호작용", **m})
print("\n=== 1) LR + 상호작용항 ===")
print(pd.Series(m).round(4))
coef1 = pd.Series(lr_inter.coef_, index=X_train.columns).sort_values(key=np.abs, ascending=False)
print(coef1.round(4))
coef1.to_csv(MODEL_DIR / "lr_interaction_coefficients.csv", header=["coefficient"])

# ---------- 2) protein_id 기준 그룹 분할 ----------
gss = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=42)
idx_train, idx_test = next(gss.split(df, groups=df["protein_id"]))
gtrain, gtest = df.iloc[idx_train], df.iloc[idx_test]
overlap = set(gtrain["protein_id"]) & set(gtest["protein_id"])
print(f"\n=== 2) GroupShuffleSplit(protein_id) === train {len(gtrain)} / test {len(gtest)}, 겹치는 protein_id 수: {len(overlap)}")

lr_group = LinearRegression().fit(gtrain[FEATURES], gtrain[TARGET])
pred_g = lr_group.predict(gtest[FEATURES])
m_g = metrics(gtest[TARGET], pred_g)
results.append({"experiment": "2) LR + protein_id 그룹 분할", **m_g})
print(pd.Series(m_g).round(4))

# ---------- 3) Ridge (StandardScaler) : 7피처 vs 9피처 ----------
FEATURES_9 = FEATURES + ["logp_pi_interaction", "mw_ratio"]

for label, feats in [("3a) Ridge 7피처", FEATURES), ("3b) Ridge 9피처(+logp_pi_interaction,+mw_ratio)", FEATURES_9)]:
    scaler = StandardScaler().fit(train_df[feats])
    Xtr = scaler.transform(train_df[feats])
    Xte = scaler.transform(test_df[feats])
    ridge = Ridge(alpha=1.0, random_state=42).fit(Xtr, train_df[TARGET])
    pred_r = ridge.predict(Xte)
    m_r = metrics(test_df[TARGET], pred_r)
    results.append({"experiment": label, **m_r})
    print(f"\n=== {label} ===")
    print(pd.Series(m_r).round(4))
    coef = pd.Series(ridge.coef_, index=feats).sort_values(key=np.abs, ascending=False)
    print("표준화 계수:")
    print(coef.round(4))
    fname = "ridge7" if "7피처" in label else "ridge9"
    coef.to_csv(MODEL_DIR / f"{fname}_coefficients_standardized.csv", header=["coefficient"])
    joblib.dump({"model": ridge, "scaler": scaler, "feature_columns": feats}, MODEL_DIR / f"{fname}.pkl")

# ---------- 요약 ----------
summary = pd.DataFrame(results)
baseline = pd.read_csv(MODEL_DIR / "linear_regression_metrics.csv", index_col=0).loc["test"]
summary_full = pd.concat([
    pd.DataFrame([{"experiment": "0) 기존 Linear Regression (baseline, 7피처)", **baseline.to_dict()}]),
    summary,
], ignore_index=True)
print("\n=== 전체 요약 (test 기준) ===")
print(summary_full.round(4).to_string(index=False))
summary_full.to_csv(MODEL_DIR / "next_steps_summary.csv", index=False)
print("\n저장 완료: next_steps_summary.csv, lr_interaction_coefficients.csv, ridge7/9_coefficients_standardized.csv, ridge7/9.pkl")
