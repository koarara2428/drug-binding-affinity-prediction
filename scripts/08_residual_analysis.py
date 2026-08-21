"""Linear Regression residual analysis — test set 기준. 이상치는 삭제하지 않고 있는 그대로 진단만 한다."""
import pandas as pd
import numpy as np
import joblib
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
SPLIT_DIR = BASE_DIR / "data" / "processed" / "split"
FIG_DIR = BASE_DIR / "output" / "day6" / "figures"
MODEL_DIR = BASE_DIR / "output" / "day6" / "models"

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

bundle = joblib.load(MODEL_DIR / "linear_regression.pkl")
model = bundle["model"]

test_df = pd.read_csv(SPLIT_DIR / "test.csv").reset_index(drop=True)
X_test, y_test = test_df[FEATURES], test_df[TARGET]
pred = model.predict(X_test)
resid = y_test - pred  # 실제값 - 예측값. 양수면 과소예측, 음수면 과대예측

res_df = test_df.copy()
res_df["predicted"] = pred
res_df["residual"] = resid
res_df["abs_residual"] = resid.abs()

# --- 1) residual distribution ---
print("=== 1) residual 분포 ===")
print(resid.describe())
print("skew:", resid.skew())

fig, ax = plt.subplots(figsize=(6.5, 4.5))
sns.histplot(resid, bins=35, color=CAT[0], ax=ax)
ax.axvline(0, color=CAT[1], linestyle="--", linewidth=1.2, label="0 (오차 없음)")
ax.axvline(resid.mean(), color=CAT[2], linestyle="--", linewidth=1.2, label=f"평균={resid.mean():.3f}")
ax.set_title("Residual 분포 (실제값 - 예측값)")
ax.set_xlabel("residual")
ax.legend()
fig.savefig(FIG_DIR / "residual_distribution.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# --- 2) predicted vs residual plot ---
fig, ax = plt.subplots(figsize=(7, 5))
ax.scatter(pred, resid, s=14, alpha=0.4, color=CAT[0])
ax.axhline(0, color=CAT[1], linestyle="--", linewidth=1.2)
ax.set_xlabel("예측값 (predicted)")
ax.set_ylabel("residual (실제값 - 예측값)")
ax.set_title("예측값 vs Residual")
fig.savefig(FIG_DIR / "predicted_vs_residual.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# --- 3) 가장 큰 오차 10개 ---
top10 = res_df.reindex(res_df["abs_residual"].sort_values(ascending=False).index).head(10)
print("\n=== 3) 오차가 가장 큰 10개 (실제값, 예측값, residual 포함) ===")
cols_order = [TARGET, "predicted", "residual"] + FEATURES
print(top10[cols_order].to_string(index=False))
top10[cols_order].to_csv(MODEL_DIR / "linear_regression_top10_errors.csv", index=False)

# --- 5) 특정 feature 범위에서 오차가 커지는지 ---
print("\n=== 5) |residual|과 각 feature의 상관계수 ===")
corr_with_abs_resid = res_df[FEATURES].apply(lambda col: col.corr(res_df["abs_residual"]))
print(corr_with_abs_resid.sort_values(key=np.abs, ascending=False))

print("\n=== feature 4분위 구간별 평균 |residual| (상위 2개 feature) ===")
top_feats = corr_with_abs_resid.abs().sort_values(ascending=False).index[:2].tolist()
for feat in top_feats:
    res_df[f"{feat}_bin"] = pd.qcut(res_df[feat], 4, duplicates="drop")
    print(f"\n--- {feat} ---")
    print(res_df.groupby(f"{feat}_bin", observed=True)["abs_residual"].agg(["mean", "count"]))

# --- 상호작용 가설 검증: logp x protein_pi 값이 클수록 오차가 커지는가 ---
res_df["logp_pi_proxy"] = res_df["logp"] * res_df["protein_pi"]
print("\n=== residual과 target(실제값)의 상관관계 (과소/과대예측 편향 확인) ===")
print("corr(residual, 실제값):", round(resid.corr(y_test), 4))
print("=== |residual|과 logp*protein_pi(제외했던 상호작용항) 상관계수 ===")
print("corr(|residual|, logp*protein_pi):", round(res_df["logp_pi_proxy"].corr(res_df["abs_residual"]), 4))

print("\n저장 완료:")
print("-", FIG_DIR / "residual_distribution.png")
print("-", FIG_DIR / "predicted_vs_residual.png")
print("-", MODEL_DIR / "linear_regression_top10_errors.csv")
