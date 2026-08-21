"""3개 모델(Linear Regression / RandomForest / XGBoost) 성능 비교"""
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
MODEL_DIR = BASE_DIR / "output" / "day6" / "models"
FIG_DIR = BASE_DIR / "output" / "day6" / "figures"

CAT = ["#2a78d6", "#eb6834", "#1baf7a"]
plt.rcParams.update({
    "font.family": "Malgun Gothic", "axes.unicode_minus": False,
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
    "axes.edgecolor": "#c3c2b7", "axes.labelcolor": "#0b0b0b", "text.color": "#0b0b0b",
    "xtick.color": "#898781", "ytick.color": "#898781",
    "axes.grid": True, "grid.color": "#e1e0d9", "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
})

files = {
    "Linear Regression": "linear_regression_metrics.csv",
    "RandomForest": "random_forest_metrics.csv",
    "XGBoost": "xgboost_metrics.csv",
}

rows = []
for name, fname in files.items():
    df = pd.read_csv(MODEL_DIR / fname, index_col=0)
    rows.append({"model": name, "split": "train", **df.loc["train"].to_dict()})
    rows.append({"model": name, "split": "test", **df.loc["test"].to_dict()})

comp = pd.DataFrame(rows)
comp["gap_R2"] = comp.groupby("model")["R2"].transform(lambda s: s.iloc[0] - s.iloc[1])
print(comp.round(4).to_string(index=False))
comp.to_csv(MODEL_DIR / "model_comparison.csv", index=False)

test_comp = comp[comp["split"] == "test"].set_index("model")
fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
for ax, metric in zip(axes, ["R2", "MAE", "RMSE"]):
    ax.bar(test_comp.index, test_comp[metric], color=CAT[:len(test_comp)])
    ax.set_title(f"test {metric}")
    for i, v in enumerate(test_comp[metric]):
        ax.text(i, v, f"{v:.3f}", ha="center", va="bottom", fontsize=9)
fig.suptitle("모델별 test 성능 비교")
fig.tight_layout()
fig.savefig(FIG_DIR / "model_comparison.png", dpi=150, bbox_inches="tight")
plt.close(fig)

print("\n저장 완료: model_comparison.csv, model_comparison.png")
