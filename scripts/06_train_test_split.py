"""STAGE ④ 4-2 — Linear Regression용 train/test 분할
피처 7개(logp, protein_pi, molecular_weight, protein_length, hydrophobicity, rotatable_bonds, polar_surface_area)
+ target(binding_affinity)만 사용. 무작위 8:2, random_state=42 (SOP 기본값, 그룹 분할 안 함).
"""
import pandas as pd
from pathlib import Path
from sklearn.model_selection import train_test_split

BASE_DIR = Path(__file__).resolve().parent.parent
DATA = BASE_DIR / "data" / "processed" / "drug_discovery_virtual_screening_processed.csv"
OUT_DIR = BASE_DIR / "data" / "processed" / "split"
OUT_DIR.mkdir(parents=True, exist_ok=True)

FEATURES = [
    "logp", "protein_pi", "molecular_weight",
    "protein_length", "hydrophobicity", "rotatable_bonds", "polar_surface_area",
]
TARGET = "binding_affinity"

df = pd.read_csv(DATA)
df = df[FEATURES + [TARGET]]

train_df, test_df = train_test_split(df, test_size=0.2, random_state=42)

print(f"전체: {len(df)}행")
print(f"train: {len(train_df)}행 ({len(train_df)/len(df)*100:.1f}%)")
print(f"test:  {len(test_df)}행 ({len(test_df)/len(df)*100:.1f}%)")

train_df.to_csv(OUT_DIR / "train.csv", index=False)
test_df.to_csv(OUT_DIR / "test.csv", index=False)

# 검증
check_train = pd.read_csv(OUT_DIR / "train.csv")
check_test = pd.read_csv(OUT_DIR / "test.csv")
assert len(check_train) + len(check_test) == len(df)
assert set(check_train.columns) == set(FEATURES + [TARGET])
assert set(check_test.columns) == set(FEATURES + [TARGET])
print("\n저장 완료:", OUT_DIR / "train.csv", "/", OUT_DIR / "test.csv")

print("\ntrain target 기술통계:")
print(train_df[TARGET].describe())
print("\ntest target 기술통계:")
print(test_df[TARGET].describe())
