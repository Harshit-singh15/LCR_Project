import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

Path(
    "Fig4/03_plots"
).mkdir(
    parents=True,
    exist_ok=True
)

df = pd.read_csv(
    "Fig4\\02_metrices\\jaccard_matrix.tsv",
    sep="\t",
    index_col=0
)

# remove suffix for cleaner labels

df.index = [
    x.replace("_sorted.bed","")
    for x in df.index
]

df.columns = [
    x.replace("_sorted.bed","")
    for x in df.columns
]

# mask upper triangle

mask = np.triu(
    np.ones_like(
        df,
        dtype=bool
    ),
    k=1
)

plt.figure(
    figsize=(10,8)
)

sns.heatmap(

    df,

    mask=mask,

    cmap="Reds",

    vmin=0,

    vmax=1,

    annot=True,

    fmt=".2f",

    square=True,

    linewidths=0.5,

    cbar_kws={
        "label":
        "Jaccard similarity"
    }
)

plt.title(
    "Pairwise overlap among LCR detection methods"
)

plt.tight_layout()

plt.savefig(
    "Fig4/03_plots/Fig4_jaccard_heatmap.png",
    dpi=600,
    bbox_inches="tight"
)

plt.show()