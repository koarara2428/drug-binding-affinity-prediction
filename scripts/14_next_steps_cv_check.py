"""13번 실험 결과(특히 Ridge 9피처의 R² 0.553 상승)가 단일 분할의 우연이 아닌지 5-fold 교차검증으로 재확인."""
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.model_selection import KFold, cross_validate
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

BASE_DIR = Path(__file__).resolve().parent.parent
PROCESSED = BASE_DIR / "data" / "processed" / "drug_discovery_virtual_screening_processed.csv"
MODEL_DIR = BASE_DIR / "output" / "day6" / "models"

FEATURES = [
    "logp", "protein_pi", "molecular_weight",
    "protein_length", "hydrophobicity", "rotatable_bonds", "polar_surface_area",
]
FEATURES_INTER = FEATURES + ["logp_x_protein_pi"]
FEATURES_9 = FEATURES + ["logp_pi_interaction", "mw_ratio"]
TARGET = "binding_affinity"

df = pd.read_csv(PROCESSED)
df["logp_x_protein_pi"] = df["logp"] * df["protein_pi"]

cv = KFold(n_splits=5, shuffle=True, random_state=42)
scoring = ["r2", "neg_mean_absolute_error", "neg_root_mean_squared_error"]

configs = {
    "0) LR 7피처 (baseline)": (LinearRegression(), FEATURES),
    "1) LR + 상호작용항": (LinearRegression(), FEATURES_INTER),
    "3a) Ridge 7피처": (make_pipeline(StandardScaler(), Ridge(alpha=1.0, random_state=42)), FEATURES),
    "3b) Ridge 9피처": (make_pipeline(StandardScaler(), Ridge(alpha=1.0, random_state=42)), FEATURES_9),
}

rows = []
for name, (model, feats) in configs.items():
    result = cross_validate(model, df[feats], df[TARGET], cv=cv, scoring=scoring)
    r2 = result["test_r2"]
    mae = -result["test_neg_mean_absolute_error"]
    rmse = -result["test_neg_root_mean_squared_error"]
    rows.append({
        "experiment": name,
        "R2_mean": r2.mean(), "R2_std": r2.std(),
        "MAE_mean": mae.mean(), "MAE_std": mae.std(),
        "RMSE_mean": rmse.mean(), "RMSE_std": rmse.std(),
    })

cv_df = pd.DataFrame(rows)
print("=== 5-fold 교차검증 (mean ± std) ===")
print(cv_df.round(4).to_string(index=False))
cv_df.to_csv(MODEL_DIR / "next_steps_cv_results.csv", index=False)

# 기존 baseline 3모델 CV 결과와 비교
prev = pd.read_csv(MODEL_DIR / "cross_validation_results.csv")
best_prev = prev.loc[prev["R2_mean"].idxmax()]
ridge9 = cv_df[cv_df["experiment"] == "3b) Ridge 9피처"].iloc[0]
diff = ridge9["R2_mean"] - best_prev["R2_mean"]
ref_std = max(ridge9["R2_std"], best_prev["R2_std"])
verdict = "유의미한 차이" if diff > ref_std else "사실상 동률"
print(f"\nRidge 9피처({ridge9['R2_mean']:.4f}) vs 기존 최고 {best_prev['model']}({best_prev['R2_mean']:.4f}): "
      f"차이={diff:.4f}, 기준std={ref_std:.4f} -> {verdict}")
