import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

INPUT_DIR = Path("Fig5/02_metrics")
OUTPUT_DIR = Path("Fig5/03_plots")

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

files = sorted(
    INPUT_DIR.glob("*_metrics.tsv")
)

n = len(files)

cols = 4
rows = int(np.ceil(n / cols))

fig, axes = plt.subplots(
    rows,
    cols,
    figsize=(16, 12)
)

axes = axes.flatten()

for ax, file in zip(axes, files):

    method = (
        file.stem
        .replace("_metrics", "")
    )

    df = pd.read_csv(
        file,
        sep="\t"
    )

    # 10 x 10 bins
    heatmap, xedges, yedges = np.histogram2d(
        df["Mutation_Percent"],
        df["Most_Frequent_AA_Percent"],
        bins=10,
        range=[
            [0,100],
            [0,100]
        ]
    )

    heatmap = np.log10(
        heatmap + 1
    )

    im = ax.imshow(
        heatmap.T,
        origin="lower",
        aspect="auto",
        extent=[
            0,100,
            0,100
        ]
    )

    ax.set_title(
        method,
        fontsize=11,
        fontweight="bold"
    )

    ax.set_xticks(
        np.arange(
            10,
            101,
            10
        )
    )

    ax.set_yticks(
        np.arange(
            10,
            101,
            10
        )
    )

    cbar = plt.colorbar(
        im,
        ax=ax,
        fraction=0.046,
        pad=0.04
    )

    cbar.set_label(
        "Count"
    )

# remove unused panels

for i in range(
    len(files),
    len(axes)
):
    fig.delaxes(
        axes[i]
    )

fig.supxlabel(
    "Mutation Percent",
    fontsize=18
)

fig.supylabel(
    "Most Frequent AA Percent",
    fontsize=18
)

plt.tight_layout()

plt.savefig(
    OUTPUT_DIR /
    "Fig5_mutation_vs_composition.png",
    dpi=600,
    bbox_inches="tight"
)

plt.show()