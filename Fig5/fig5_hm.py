import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pathlib import Path

# ======================================================
# Input / Output
# ======================================================

INPUT_DIR = Path("celegans\\dataforFig5\\02_metrics")
OUTPUT_DIR = Path("celegans\\Fig_outputs\\Fig5")

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

# ======================================================
# Clean method name
# ======================================================

organisms = [
    "_mouse",
    "_celegans",
    "_arabidopsis",
    "_ecoli"
]

def clean_name(name):

    name = name.replace("_metrics", "")

    for org in organisms:
        name = name.replace(org, "")

    return name

# ======================================================
# Plot every tool
# ======================================================

for file in sorted(INPUT_DIR.glob("*_metrics.tsv")):

    method = clean_name(file.stem)

    print(f"Processing {method}")

    df = pd.read_csv(
        file,
        sep="\t"
    )

    # ------------------------------------------
    # Bin values exactly like paper
    # ------------------------------------------

    df["MutBin"] = (
        df["Mutation_Percent"] / 10
    ).round().astype(int)

    df["FreqBin"] = (
        df["Most_Common_AA_Percent"] / 10
    ).round().astype(int)

    # Keep within plotting range

    df["MutBin"] = df["MutBin"].clip(0,10)

    df["FreqBin"] = df["FreqBin"].clip(0,10)

    # ------------------------------------------
    # Count frequencies
    # ------------------------------------------

    heatmap = np.zeros((11,11))

    for _, row in df.iterrows():

        heatmap[
            row["FreqBin"],
            row["MutBin"]
        ] += 1

    # ------------------------------------------
    # Plot
    # ------------------------------------------

    fig, ax = plt.subplots(
        figsize=(7,6)
    )

    vmax = max(
        1,
        heatmap.max()
    )

    im = ax.imshow(

        heatmap,

        origin="lower",

        cmap="Reds",

        norm=LogNorm(
            vmin=1,
            vmax=vmax
        ),

        aspect="auto"

    )

    # ------------------------------------------

    ax.set_xticks(range(11))
    ax.set_yticks(range(11))

    ax.set_xticklabels(
        [f"{i*10}" for i in range(11)],
        rotation=45
    )

    ax.set_yticklabels(
        [f"{i*10}" for i in range(11)]
    )

    ax.set_xlabel(
        "Mutation Percentage"
    )

    ax.set_ylabel(
        "Most Frequent AA Percentage"
    )

    ax.set_title(
        method
    )

    cbar = plt.colorbar(
        im,
        ax=ax
    )

    cbar.set_label(
        "Number of LCRs"
    )

    plt.tight_layout()

    # ------------------------------------------
    # Save
    # ------------------------------------------

    plt.savefig(
        OUTPUT_DIR /
        f"{method}_Fig5.png",
        dpi=600,
        bbox_inches="tight"
    )


    plt.close()

print("Done.")