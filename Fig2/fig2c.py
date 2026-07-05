from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

df = pd.read_csv(
    Path(sys.argv[1]),
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
    "magma",
    n_colors=13
)

# FIX 1: Assigned 'x' ("Consensus") to 'hue' and set 'legend=False'
ax = sns.boxplot(
    data=df,
    x="Consensus",
    y="Purity",
    hue="Consensus",
    palette=palette,
    legend=False,
    showfliers=False,
    linewidth=1
)

plt.xlabel("Consensus level")
plt.ylabel("Purity")

# FIX 2: Explicitly set the tick positions before setting the tick labels
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels)

plt.title(
    "Fig 2C : Compositional purity across consensus levels"
)

plt.tight_layout()

plt.savefig(
    Path(sys.argv[2]) / "Fig2C_purity.png",
    dpi=600,
    bbox_inches="tight"
)
