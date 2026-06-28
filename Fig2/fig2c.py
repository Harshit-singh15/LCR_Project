import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv(
    r"fruitfly\dataforFig2\purity\purity.tsv",
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
    "Fig 2C: Fruit Fly : Compositional purity across consensus levels"
)

plt.tight_layout()

plt.savefig(
    r"fruitfly\Fig_outputs\Fig2\Fig2C_purity.png",
    dpi=600,
    bbox_inches="tight"
)

plt.show()