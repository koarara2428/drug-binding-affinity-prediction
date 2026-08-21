"""STAGE ① 원본 탐색 + STAGE ③ 핵심 내용 시각화(일부) — data/raw/drug_discovery_virtual_screening.csv
읽기 전용. data/raw/는 수정하지 않는다.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
RAW = BASE_DIR / "data" / "raw" / "drug_discovery_virtual_screening.csv"
FIG_DIR = BASE_DIR / "output" / "day6" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

CAT = ["#2a78d6", "#eb6834", "#1baf7a"]
DIVERGING = LinearSegmentedColormap.from_list("div", ["#2a78d6", "#f0efec", "#e34948"])
plt.rcParams.update({
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
    "axes.edgecolor": "#c3c2b7", "axes.labelcolor": "#0b0b0b", "text.color": "#0b0b0b",
    "xtick.color": "#898781", "ytick.color": "#898781",
    "axes.grid": True, "grid.color": "#e1e0d9", "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
})
sns.set_palette(CAT)

df = pd.read_csv(RAW)
numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()

print("=== shape ===")
print(df.shape)

print("\n=== dtypes ===")
print(df.dtypes)

print("\n=== 결측 비율(%) ===")
print((df.isna().mean() * 100).round(2))

print("\n=== describe (수치형) ===")
print(df[numeric_cols].describe().T.round(3))

# --- 이상치: IQR 기준 (그룹 변수 protein_id는 400개라 그룹별 IQR은 표본 부족으로 생략) ---
print("\n=== IQR 이상치 개수 (전체 기준) ===")
for col in numeric_cols:
    q1, q3 = df[col].quantile([0.25, 0.75])
    iqr = q3 - q1
    lo, hi = q1 - 1.5 * iqr, q3 + 1.5 * iqr
    n_out = ((df[col] < lo) | (df[col] > hi)).sum()
    if n_out > 0:
        print(f"{col}: {n_out}건 (범위 [{lo:.2f}, {hi:.2f}])")

# --- 다중공선성 ---
corr = df[numeric_cols].corr()
c_abs = corr.abs()
pairs = c_abs.where(~np.tril(np.ones(c_abs.shape)).astype(bool)).stack().sort_values(ascending=False)
print("\n=== 다중공선성 (|r| >= 0.8) ===")
print(pairs[pairs >= 0.8])

# --- 선형성 참고용 상관계수 (binding_affinity vs 주요 피처) ---
print("\n=== binding_affinity와의 상관계수 ===")
print(corr["binding_affinity"].sort_values(ascending=False))

# --- 시각화 1: 상관 히트맵 ---
fig, ax = plt.subplots(figsize=(9, 7))
sns.heatmap(corr, annot=True, fmt=".2f", cmap=DIVERGING, center=0, vmin=-1, vmax=1,
            linewidths=2, linecolor="#fcfcfb", ax=ax, annot_kws={"size": 7})
ax.set_title("Correlation Heatmap (numeric columns)")
fig.savefig(FIG_DIR / "correlation_heatmap.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# --- 시각화 2: binding_affinity 분포 ---
fig, ax = plt.subplots(figsize=(6, 4))
sns.histplot(df["binding_affinity"], bins=40, color=CAT[0], ax=ax)
ax.axvline(7.0, color=CAT[1], linestyle="--", linewidth=1, label="active threshold (7.0)")
ax.set_title("binding_affinity distribution")
ax.legend()
fig.savefig(FIG_DIR / "histogram_binding_affinity.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# --- 시각화 3: target vs 주요 피처 산점도(선형성 확인) ---
features = ["molecular_weight", "logp", "protein_pi"]
fig, axes = plt.subplots(1, len(features), figsize=(4.5 * len(features), 4))
for ax, feat in zip(axes, features):
    ax.scatter(df[feat], df["binding_affinity"], s=8, alpha=0.35, color=CAT[0])
    r = df[feat].corr(df["binding_affinity"])
    ax.set_title(f"{feat} vs binding_affinity (r={r:.2f})")
    ax.set_xlabel(feat)
    ax.set_ylabel("binding_affinity")
fig.tight_layout()
fig.savefig(FIG_DIR / "scatter_target_vs_features.png", dpi=150, bbox_inches="tight")
plt.close(fig)

print("\n저장된 그림:", [p.name for p in FIG_DIR.glob("*.png")])
