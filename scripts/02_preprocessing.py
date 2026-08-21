"""STAGE ② 전처리 — data/raw/drug_discovery_virtual_screening.csv
- active, compound_id 제외 (누수/식별자)
- logp, polar_surface_area, hydrophobicity 결측 행 삭제
원본(data/raw/)은 수정하지 않는다. 결과는 data/processed/에 저장.
"""
import pandas as pd
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
RAW = BASE_DIR / "data" / "raw" / "drug_discovery_virtual_screening.csv"
OUT = BASE_DIR / "data" / "processed" / "drug_discovery_virtual_screening_processed.csv"
OUT.parent.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(RAW)
n_before = len(df)

DROP_COLS = ["active", "compound_id"]
MISSING_SUBSET = ["logp", "polar_surface_area", "hydrophobicity"]

df = df.drop(columns=DROP_COLS)

# polar_surface_area는 정의상 항상 0 이상. 음수는 물리적으로 불가능한 값이라 결측으로 재분류 후 삭제한다.
n_impossible = (df["polar_surface_area"] < 0).sum()
df.loc[df["polar_surface_area"] < 0, "polar_surface_area"] = pd.NA
print(f"polar_surface_area 음수(물리적으로 불가능) {n_impossible}건을 결측으로 재분류")

df = df.dropna(subset=MISSING_SUBSET).reset_index(drop=True)

n_after = len(df)
print(f"행 수: {n_before} -> {n_after} ({n_before - n_after}행 삭제, {(n_before-n_after)/n_before*100:.1f}%)")
print(f"컬럼 수: {df.shape[1]} (제외: {DROP_COLS})")

df.to_csv(OUT, index=False)

check = pd.read_csv(OUT)
assert check.shape == df.shape, "저장 후 shape 불일치"
assert check["protein_id"].notna().all()
assert check[MISSING_SUBSET].isna().sum().sum() == 0, "결측이 남아있음"
assert (check["polar_surface_area"] >= 0).all(), "음수 polar_surface_area가 남아있음"
assert not set(DROP_COLS) & set(check.columns), "제외 컬럼이 남아있음"
print("검증 통과:", OUT)
print(check.dtypes)
