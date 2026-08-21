"""주요 수치형 feature 전체 분포 확인 — data/processed 기준. binding_affinity는 이미 별도 EDA 완료라 제외."""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
DATA = BASE_DIR / "data" / "processed" / "drug_discovery_virtual_screening_processed.csv"
FIG_DIR = BASE_DIR / "output" / "day6" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

CAT = ["#2a78d6", "#eb6834", "#1baf7a"]
plt.rcParams.update({
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
    "axes.edgecolor": "#c3c2b7", "axes.labelcolor": "#0b0b0b", "text.color": "#0b0b0b",
    "xtick.color": "#898781", "ytick.color": "#898781",
    "axes.grid": True, "grid.color": "#e1e0d9", "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
})

df = pd.read_csv(DATA)
cols = [c for c in df.select_dtypes(include=[np.number]).columns if c != "binding_affinity"]

# --- 히스토그램 그리드 (13개 feature) ---
n = len(cols)
ncols = 4
nrows = -(-n // ncols)
fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 3.2 * nrows))
axes = axes.flatten()
for ax, col in zip(axes, cols):
    sns.histplot(df[col], bins=35, color=CAT[0], ax=ax)
    ax.set_title(col, fontsize=10)
    ax.set_xlabel("")
for ax in axes[n:]:
    ax.axis("off")
fig.tight_layout()
fig.savefig(FIG_DIR / "numeric_features_histograms.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# --- boxplot: 이상값 비율 상위 4개(outlier% > 1%) ---
notable = ["mw_ratio", "rotatable_bonds", "logp", "molecular_weight"]
fig, axes = plt.subplots(1, len(notable), figsize=(3.6 * len(notable), 4.2))
for ax, col in zip(axes, notable):
    sns.boxplot(y=df[col], color=CAT[0], ax=ax, width=0.35)
    ax.set_title(col, fontsize=10)
fig.tight_layout()
fig.savefig(FIG_DIR / "numeric_features_boxplots_notable.png", dpi=150, bbox_inches="tight")
plt.close(fig)

print("저장 완료:", ["numeric_features_histograms.png", "numeric_features_boxplots_notable.png"])
