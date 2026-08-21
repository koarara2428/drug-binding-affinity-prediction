"""README Next steps 4) 최종 모델 재학습
CV로 확인된 최고 성능 구성(Ridge, 9피처: 7개 원본 + logp_pi_interaction + mw_ratio, StandardScaler)을
train+test 전체(1,825행)로 재학습해 최종 모델로 저장.
"""
import pandas as pd
import numpy as np
import joblib
from pathlib import Path
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error

BASE_DIR = Path(__file__).resolve().parent.parent
PROCESSED = BASE_DIR / "data" / "processed" / "drug_discovery_virtual_screening_processed.csv"
MODEL_DIR = BASE_DIR / "output" / "day6" / "models"

FEATURES_FINAL = [
    "logp", "protein_pi", "molecular_weight", "protein_length",
    "hydrophobicity", "rotatable_bonds", "polar_surface_area",
    "logp_pi_interaction", "mw_ratio",
]
TARGET = "binding_affinity"

df = pd.read_csv(PROCESSED)
X, y = df[FEATURES_FINAL], df[TARGET]

final_model = make_pipeline(StandardScaler(), Ridge(alpha=1.0, random_state=42))
final_model.fit(X, y)

pred = final_model.predict(X)
m = {
    "R2": r2_score(y, pred),
    "MAE": mean_absolute_error(y, pred),
    "RMSE": np.sqrt(mean_squared_error(y, pred)),
}
print("=== 전체 데이터(1,825행) 재학습 후 자체 적합도 (참고용, 일반화 성능 아님) ===")
print(pd.Series(m).round(4))
print("주의: 이 수치는 학습에 쓴 데이터 자체를 다시 맞춘 것이라 과적합 여부를 알 수 없음.")
print("일반화 성능은 앞서 5-fold CV 결과(R2 0.571 ± 0.068)를 참고할 것.")

ridge_step = final_model.named_steps["ridge"]
coef = pd.Series(ridge_step.coef_, index=FEATURES_FINAL).sort_values(key=np.abs, ascending=False)
print("\n표준화 계수:")
print(coef.round(4))
coef.to_csv(MODEL_DIR / "final_model_coefficients.csv", header=["coefficient"])

joblib.dump({"model": final_model, "feature_columns": FEATURES_FINAL, "cv_r2_mean": 0.5714, "cv_r2_std": 0.0681},
            MODEL_DIR / "final_model.pkl")
print("\n저장 완료:", MODEL_DIR / "final_model.pkl", "/", MODEL_DIR / "final_model_coefficients.csv")
