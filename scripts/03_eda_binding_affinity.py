"""binding_affinity EDA — data/processed/drug_discovery_virtual_screening_processed.csv 기준"""
import pandas as pd
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
s = df["binding_affinity"]

stats = {
    "n": len(s), "mean": s.mean(), "median": s.median(), "std": s.std(),
    "min": s.min(), "max": s.max(), "skew": s.skew(),
}
for k, v in stats.items():
    print(f"{k}: {v}")

q1, q3 = s.quantile([0.25, 0.75])
iqr = q3 - q1
lo, hi = q1 - 1.5 * iqr, q3 + 1.5 * iqr

fig, ax = plt.subplots(figsize=(7, 4.5))
sns.histplot(s, bins=40, color=CAT[0], ax=ax)
ax.axvline(stats["mean"], color=CAT[1], linestyle="--", linewidth=1.2, label=f"mean={stats['mean']:.2f}")
ax.axvline(stats["median"], color=CAT[2], linestyle="--", linewidth=1.2, label=f"median={stats['median']:.2f}")
ax.set_title("binding_affinity histogram (processed, n=1826)")
ax.legend()
fig.savefig(FIG_DIR / "binding_affinity_histogram.png", dpi=150, bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(6, 4.5))
sns.boxplot(y=s, color=CAT[0], ax=ax, width=0.3)
ax.set_title("binding_affinity boxplot (processed, n=1826)")
ax.set_ylabel("binding_affinity")
fig.savefig(FIG_DIR / "binding_affinity_boxplot.png", dpi=150, bbox_inches="tight")
plt.close(fig)

print("정상범위(IQR 1.5배):", lo, hi)
print("이상값 후보:", ((s < lo) | (s > hi)).sum())
