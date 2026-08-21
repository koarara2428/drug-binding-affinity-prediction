"""STAGE ④ 4-7 — 교차검증 (5-fold). train+test 전체(1,825행)를 다시 합쳐서 사용.
목적: RandomForest vs XGBoost의 test 성능 차이(0.477 vs 0.472)가 진짜 차이인지, 단일 분할의 우연인지 확인.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.model_selection import KFold, cross_validate
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor

BASE_DIR = Path(__file__).resolve().parent.parent
SPLIT_DIR = BASE_DIR / "data" / "processed" / "split"
MODEL_DIR = BASE_DIR / "output" / "day6" / "models"
FIG_DIR = BASE_DIR / "output" / "day6" / "figures"

CAT = ["#2a78d6", "#eb6834", "#1baf7a"]
plt.rcParams.update({
    "font.family": "Malgun Gothic", "axes.unicode_minus": False,
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
    "axes.edgecolor": "#c3c2b7", "axes.labelcolor": "#0b0b0b", "text.color": "#0b0b0b",
    "xtick.color": "#898781", "ytick.color": "#898781",
    "axes.grid": True, "grid.color": "#e1e0d9", "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
})

FEATURES = [
    "logp", "protein_pi", "molecular_weight",
    "protein_length", "hydrophobicity", "rotatable_bonds", "polar_surface_area",
]
TARGET = "binding_affinity"

train_df = pd.read_csv(SPLIT_DIR / "train.csv")
test_df = pd.read_csv(SPLIT_DIR / "test.csv")
df_all = pd.concat([train_df, test_df], ignore_index=True)
X, y = df_all[FEATURES], df_all[TARGET]
print(f"전체 표본: {len(df_all)}행 (train {len(train_df)} + test {len(test_df)})")

cv = KFold(n_splits=5, shuffle=True, random_state=42)
scoring = ["r2", "neg_mean_absolute_error", "neg_root_mean_squared_error"]

models = {
    "Linear Regression": LinearRegression(),
    "RandomForest": RandomForestRegressor(n_estimators=300, max_depth=6, min_samples_leaf=8, random_state=42, n_jobs=-1),
    "XGBoost": XGBRegressor(n_estimators=300, max_depth=3, learning_rate=0.05, reg_alpha=0.1, reg_lambda=1.0,
                             subsample=0.8, colsample_bytree=0.8, random_state=42, n_jobs=-1),
}

rows = []
raw_scores = {}
for name, model in models.items():
    result = cross_validate(model, X, y, cv=cv, scoring=scoring)
    r2 = result["test_r2"]
    mae = -result["test_neg_mean_absolute_error"]
    rmse = -result["test_neg_root_mean_squared_error"]
    raw_scores[name] = r2
    rows.append({
        "model": name,
        "R2_mean": r2.mean(), "R2_std": r2.std(),
        "MAE_mean": mae.mean(), "MAE_std": mae.std(),
        "RMSE_mean": rmse.mean(), "RMSE_std": rmse.std(),
    })

cv_df = pd.DataFrame(rows)
print("\n=== 5-fold 교차검증 결과 (mean ± std) ===")
print(cv_df.round(4).to_string(index=False))
cv_df.to_csv(MODEL_DIR / "cross_validation_results.csv", index=False)

# 판정 규칙: mean 차이 < std -> 사실상 동률 / mean 차이 > std -> 유의미한 차이
print("\n=== 모델 쌍별 판정 (SOP 4-7 규칙) ===")
names = list(models.keys())
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        a, b = names[i], names[j]
        ra, rb = cv_df.loc[cv_df.model == a, "R2_mean"].iloc[0], cv_df.loc[cv_df.model == b, "R2_mean"].iloc[0]
        sa, sb = cv_df.loc[cv_df.model == a, "R2_std"].iloc[0], cv_df.loc[cv_df.model == b, "R2_std"].iloc[0]
        diff = abs(ra - rb)
        ref_std = max(sa, sb)
        verdict = "유의미한 차이" if diff > ref_std else "사실상 동률"
        print(f"{a}({ra:.4f}) vs {b}({rb:.4f}): 차이={diff:.4f}, 기준std={ref_std:.4f} -> {verdict}")

# boxplot: fold별 R2 분포
fig, ax = plt.subplots(figsize=(7, 5))
ax.boxplot([raw_scores[n] for n in names], tick_labels=names, patch_artist=True,
           boxprops=dict(facecolor=CAT[0], alpha=0.6))
for i, n in enumerate(names):
    ax.scatter([i + 1] * len(raw_scores[n]), raw_scores[n], color=CAT[1], zorder=3, s=25)
ax.set_ylabel("R² (5-fold)")
ax.set_title("5-fold 교차검증 R² 분포")
fig.savefig(FIG_DIR / "cross_validation_r2_boxplot.png", dpi=150, bbox_inches="tight")
plt.close(fig)

print("\n저장 완료: cross_validation_results.csv, cross_validation_r2_boxplot.png")
