import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

Path("Fig2/05_plots").mkdir(
    parents=True,
    exist_ok=True
)

df = pd.read_csv(
    "Fig2/04_metrics/entropy.tsv",
    sep="\t"
)

counts = (
    df.groupby("Consensus")
      .size()
      .to_dict()
)

labels = [
    f"{i}\n(n={counts.get(i,0)})"
    for i in range(1,14)
]

plt.figure(figsize=(12,6))

palette = sns.color_palette(
    "viridis",
    n_colors=13
)

ax = sns.boxplot(
    data=df,
    x="Consensus",
    y="Entropy",
    palette=palette,
    showfliers=False,
    linewidth=1
)

plt.xlabel("Consensus level")
plt.ylabel("Shannon entropy")

ax.set_xticklabels(labels)

plt.title(
    "Sequence complexity across consensus levels"
)

plt.tight_layout()

plt.savefig(
    "Fig2/05_plots/Fig2B_entropy.png",
    dpi=600,
    bbox_inches="tight"
)

plt.show()