"""STAGE ④ — RandomForest
target: binding_affinity / feature: Linear Regression과 동일 7개
표본이 수백 행 수준(train 1,460)이라 SOP 4-4 권장대로 규제 파라미터를 보수적으로 설정.
"""
import pandas as pd
import numpy as np
import joblib
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error

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

train_df = pd.read_csv(SPLIT_DIR / "train.csv")
test_df = pd.read_csv(SPLIT_DIR / "test.csv")
X_train, y_train = train_df[FEATURES], train_df[TARGET]
X_test, y_test = test_df[FEATURES], test_df[TARGET]

model = RandomForestRegressor(
    n_estimators=300, max_depth=6, min_samples_leaf=8,
    oob_score=True, random_state=42, n_jobs=-1,
)
model.fit(X_train, y_train)

pred_train = model.predict(X_train)
pred_test = model.predict(X_test)

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
print("\nOOB score (R2):", round(model.oob_score_, 4))
metrics_df.to_csv(MODEL_DIR / "random_forest_metrics.csv")

# 실제값 vs 예측값
fig, ax = plt.subplots(figsize=(6, 6))
ax.scatter(y_test, pred_test, s=14, alpha=0.4, color=CAT[0], label="test 예측")
lims = [min(y_test.min(), pred_test.min()) - 0.5, max(y_test.max(), pred_test.max()) + 0.5]
ax.plot(lims, lims, color=CAT[1], linestyle="--", linewidth=1.2, label="완벽한 예측선(y=x)")
ax.set_xlim(lims); ax.set_ylim(lims)
ax.set_xlabel("실제값 (binding_affinity)")
ax.set_ylabel("예측값 (binding_affinity)")
ax.set_title(f"RandomForest: 실제값 vs 예측값 (test, R²={m_test['R2']:.3f})")
ax.legend()
ax.set_aspect("equal", adjustable="box")
fig.savefig(FIG_DIR / "random_forest_actual_vs_pred.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# feature importance
imp_df = pd.DataFrame({"feature": FEATURES, "importance": model.feature_importances_}).sort_values("importance", ascending=False)
print("\n=== feature importance ===")
print(imp_df.to_string(index=False))
imp_df.to_csv(MODEL_DIR / "random_forest_feature_importance.csv", index=False)

# leaf당 평균 샘플 수 (SOP 4-6.3, 암기 여부 확인)
leaves = [t.get_n_leaves() for t in model.estimators_]
avg_samples_per_leaf = len(X_train) / (sum(leaves) / len(leaves))
print("\nleaf당 평균 샘플 수:", round(avg_samples_per_leaf, 2))

joblib.dump({"model": model, "feature_columns": FEATURES}, MODEL_DIR / "random_forest.pkl")
print("\n저장 완료: random_forest_metrics.csv, random_forest_feature_importance.csv, random_forest.pkl, random_forest_actual_vs_pred.png")
