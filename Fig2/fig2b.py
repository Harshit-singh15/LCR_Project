import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Create output directory
Path(r"fruitfly\Fig_outputs\Fig2").mkdir(
    parents=True,
    exist_ok=True
)

# Read data
df = pd.read_csv(
    r"fruitfly\dataforFig2\entropy\entropy.tsv",
    sep="\t"
)

# Get group sizes for the labels
counts = (
    df.groupby("Consensus")
      .size()
      .to_dict()
)

# Generate labels (1 to 13)
labels = [
    f"{i}\n(n={counts.get(i,0)})"
    for i in range(1, 14)
]

plt.figure(figsize=(12, 6))

palette = sns.color_palette(
    "viridis",
    n_colors=13
)

# FIXED: Added 'hue' and 'legend=False' to resolve the FutureWarning
ax = sns.boxplot(
    data=df,
    x="Consensus",
    y="Entropy",
    hue="Consensus",
    palette=palette,
    showfliers=False,
    linewidth=1,
    legend=False
)

plt.xlabel("Consensus level")
plt.ylabel("Shannon entropy")

# FIXED: Explicitly set the ticks first to resolve the UserWarning
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels)

plt.title(
    "Fruit Fly : Shannon entropy distribution across consensus levels"
)

plt.tight_layout()

plt.savefig(
    r"fruitfly\Fig_outputs\Fig2\Fig2B_entropy.png",
    dpi=600,
    bbox_inches="tight"
)

plt.show()