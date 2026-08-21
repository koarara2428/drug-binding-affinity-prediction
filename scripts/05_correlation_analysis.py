"""binding_affinity와 수치형 feature 상관관계 분석 — data/processed 최종본(1,825행) 기준"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
DATA = BASE_DIR / "data" / "processed" / "drug_discovery_virtual_screening_processed.csv"
FIG_DIR = BASE_DIR / "output" / "day6" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

DIVERGING = LinearSegmentedColormap.from_list("div", ["#2a78d6", "#f0efec", "#e34948"])
plt.rcParams.update({
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
    "axes.edgecolor": "#c3c2b7", "axes.labelcolor": "#0b0b0b", "text.color": "#0b0b0b",
    "xtick.color": "#898781", "ytick.color": "#898781",
    "axes.grid": True, "grid.color": "#e1e0d9", "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
})

df = pd.read_csv(DATA)
numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()

# 1) correlation matrix
corr = df[numeric_cols].corr()
corr.to_csv(BASE_DIR / "output" / "day6" / "correlation_matrix.csv")
print("=== correlation matrix 저장: output/day6/correlation_matrix.csv ===")

# 2) heatmap
fig, ax = plt.subplots(figsize=(9, 7.5))
sns.heatmap(corr, annot=True, fmt=".2f", cmap=DIVERGING, center=0, vmin=-1, vmax=1,
            linewidths=2, linecolor="#fcfcfb", ax=ax, annot_kws={"size": 7})
ax.set_title("Correlation Heatmap — processed (n=1825)")
fig.savefig(FIG_DIR / "correlation_heatmap_final.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# 3) binding_affinity correlation ranking
rank = corr["binding_affinity"].drop("binding_affinity").sort_values(key=np.abs, ascending=False)
print("\n=== binding_affinity와의 상관계수 순위(절대값 기준) ===")
for col, r in rank.items():
    print(f"{col:22s} r = {r:+.4f}")

fig, ax = plt.subplots(figsize=(7, 5.5))
colors = ["#e34948" if v > 0 else "#2a78d6" for v in rank.values]
ax.barh(rank.index[::-1], rank.values[::-1], color=colors[::-1])
ax.axvline(0, color="#898781", linewidth=0.8)
ax.set_title("binding_affinity correlation ranking")
ax.set_xlabel("Pearson r")
fig.savefig(FIG_DIR / "binding_affinity_correlation_ranking.png", dpi=150, bbox_inches="tight")
plt.close(fig)

print("\n저장 완료: correlation_heatmap_final.png, binding_affinity_correlation_ranking.png")
