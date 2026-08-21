"""최종 모델(Ridge 9피처) residual analysis.
주의: final_model.pkl은 train+test 전체(1,825행)로 재학습돼 held-out 데이터가 없다.
따라서 held-out 평가가 가능한 ridge9.pkl(=train만으로 학습, script 13)을 test set에 적용해 진단한다.
이상값은 삭제하지 않고 진단만 한다.
"""
import pandas as pd
import numpy as np
import joblib
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.model_selection import train_test_split

BASE_DIR = Path(__file__).resolve().parent.parent
PROCESSED = BASE_DIR / "data" / "processed" / "drug_discovery_virtual_screening_processed.csv"
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

TARGET = "binding_affinity"

bundle = joblib.load(MODEL_DIR / "ridge9.pkl")
model, scaler, FEATURES = bundle["model"], bundle["scaler"], bundle["feature_columns"]

df = pd.read_csv(PROCESSED)
_, test_df = train_test_split(df, test_size=0.2, random_state=42)  # script 13과 동일 재현
test_df = test_df.reset_index(drop=True)

# ridge9.pkl은 scaler와 model을 파이프라인이 아니라 따로 저장했으므로 직접 순서대로 적용해야 한다
X_test_scaled = scaler.transform(test_df[FEATURES])
pred = model.predict(X_test_scaled)
resid = test_df[TARGET] - pred

res_df = test_df.copy()
res_df["predicted"] = pred
res_df["residual"] = resid
res_df["abs_residual"] = resid.abs()

print("=== Ridge 9피처(held-out test, n=%d) residual 분포 ===" % len(test_df))
print(resid.describe())
print("skew:", resid.skew())
print("\ncorr(residual, 실제값):", round(resid.corr(test_df[TARGET]), 4), " (참고: 기존 LR 7피처는 0.773)")

fig, ax = plt.subplots(figsize=(6.5, 4.5))
sns.histplot(resid, bins=35, color=CAT[0], ax=ax)
ax.axvline(0, color=CAT[1], linestyle="--", linewidth=1.2, label="0 (오차 없음)")
ax.axvline(resid.mean(), color=CAT[2], linestyle="--", linewidth=1.2, label=f"평균={resid.mean():.3f}")
ax.set_title("최종 모델(Ridge 9피처) Residual 분포")
ax.set_xlabel("residual")
ax.legend()
fig.savefig(FIG_DIR / "final_model_residual_distribution.png", dpi=150, bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(7, 5))
ax.scatter(pred, resid, s=14, alpha=0.4, color=CAT[0])
ax.axhline(0, color=CAT[1], linestyle="--", linewidth=1.2)
ax.set_xlabel("예측값 (predicted)")
ax.set_ylabel("residual (실제값 - 예측값)")
ax.set_title("최종 모델(Ridge 9피처) 예측값 vs Residual")
fig.savefig(FIG_DIR / "final_model_predicted_vs_residual.png", dpi=150, bbox_inches="tight")
plt.close(fig)

top10 = res_df.reindex(res_df["abs_residual"].sort_values(ascending=False).index).head(10)
cols_order = [TARGET, "predicted", "residual"] + FEATURES
print("\n=== 오차 최대 10건 ===")
print(top10[cols_order].to_string(index=False))
top10[cols_order].to_csv(MODEL_DIR / "final_model_top10_errors.csv", index=False)

print("\n=== |residual|과 각 feature 상관계수 ===")
corr_abs = res_df[FEATURES].apply(lambda c: c.corr(res_df["abs_residual"]))
print(corr_abs.sort_values(key=np.abs, ascending=False))

# 실제값 구간별(저/중/고) 평균 residual — 극단값 편향이 남아있는지 직접 확인
res_df["target_bin"] = pd.qcut(res_df[TARGET], [0, 0.1, 0.9, 1.0], labels=["하위10%", "중간80%", "상위10%"])
print("\n=== 실제값 구간별 평균 residual (양수=과소예측, 음수=과대예측) ===")
print(res_df.groupby("target_bin", observed=True)["residual"].agg(["mean", "count"]))

print("\n저장 완료: final_model_residual_distribution.png, final_model_predicted_vs_residual.png, final_model_top10_errors.csv")
