"""STAGE ④ — Linear Regression 베이스라인
target: binding_affinity
feature: logp, protein_pi, molecular_weight, protein_length, hydrophobicity, rotatable_bonds, polar_surface_area
"""
import pandas as pd
import numpy as np
import joblib
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error

BASE_DIR = Path(__file__).resolve().parent.parent
SPLIT_DIR = BASE_DIR / "data" / "processed" / "split"
FIG_DIR = BASE_DIR / "output" / "day6" / "figures"
MODEL_DIR = BASE_DIR / "output" / "day6" / "models"
FIG_DIR.mkdir(parents=True, exist_ok=True)
MODEL_DIR.mkdir(parents=True, exist_ok=True)

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

X_train, y_train = train_df[FEATURES], train_df[TARGET]
X_test, y_test = test_df[FEATURES], test_df[TARGET]

# 1) train set으로 학습
model = LinearRegression()
model.fit(X_train, y_train)

# 2) test set으로 예측
pred_test = model.predict(X_test)
pred_train = model.predict(X_train)  # 과적합 비교용 참고치

# 3~5) R², MAE, RMSE
def metrics(y_true, y_pred):
    return {
        "R2": r2_score(y_true, y_pred),
        "MAE": mean_absolute_error(y_true, y_pred),
        "RMSE": np.sqrt(mean_squared_error(y_true, y_pred)),
    }

m_train = metrics(y_train, pred_train)
m_test = metrics(y_test, pred_test)

metrics_df = pd.DataFrame([m_train, m_test], index=["train", "test"])
print("=== 평가지표 ===")
print(metrics_df.round(4))
metrics_df.to_csv(MODEL_DIR / "linear_regression_metrics.csv")

# 6) 실제값 vs 예측값 scatter plot (test)
fig, ax = plt.subplots(figsize=(6, 6))
ax.scatter(y_test, pred_test, s=14, alpha=0.4, color=CAT[0], label="test 예측")
lims = [min(y_test.min(), pred_test.min()) - 0.5, max(y_test.max(), pred_test.max()) + 0.5]
ax.plot(lims, lims, color=CAT[1], linestyle="--", linewidth=1.2, label="완벽한 예측선(y=x)")
ax.set_xlim(lims)
ax.set_ylim(lims)
ax.set_xlabel("실제값 (binding_affinity)")
ax.set_ylabel("예측값 (binding_affinity)")
ax.set_title(f"Linear Regression: 실제값 vs 예측값 (test, R²={m_test['R2']:.3f})")
ax.legend()
ax.set_aspect("equal", adjustable="box")
fig.savefig(FIG_DIR / "linear_regression_actual_vs_pred.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# 7) coefficient 확인
coef_df = pd.DataFrame({
    "feature": FEATURES,
    "coefficient": model.coef_,
}).sort_values("coefficient", key=np.abs, ascending=False)
coef_df.loc[len(coef_df)] = ["intercept", model.intercept_]
print("\n=== 회귀계수 ===")
print(coef_df.to_string(index=False))
coef_df.to_csv(MODEL_DIR / "linear_regression_coefficients.csv", index=False)

joblib.dump({"model": model, "feature_columns": FEATURES}, MODEL_DIR / "linear_regression.pkl")
print("\n저장 완료:")
print("-", MODEL_DIR / "linear_regression_metrics.csv")
print("-", MODEL_DIR / "linear_regression_coefficients.csv")
print("-", MODEL_DIR / "linear_regression.pkl")
print("-", FIG_DIR / "linear_regression_actual_vs_pred.png")
