import re
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

# ======================================================
# Input / Output
# ======================================================

INPUT_DIR = Path(
    "celegans/dataforFig7"
)

OUTPUT_DIR = Path(
    "celegans/Fig_outputs/Fig7"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

# ======================================================
# Files
# ======================================================

panels = [

    (
        "gene_length_summary.tsv",
        "A. Gene Length"
    ),

    (
        "lcr_count_summary.tsv",
        "B. LCR Count"
    ),

    (
        "coverage_summary.tsv",
        "C. Coverage"
    ),

    (
        "entropy_ratio_summary.tsv",
        "D. Entropy Ratio"
    )

]

# ======================================================
# Plot
# ======================================================

fig, axes = plt.subplots(
    2,
    2,
    figsize=(16,12)
)

axes = axes.flatten()

legend_handles = None
legend_labels = None

for ax, (filename, title) in zip(
    axes,
    panels
):

    file = INPUT_DIR / filename

    if not file.exists():

        print(f"Skipping {filename}")

        ax.axis("off")

        continue

    df = pd.read_csv(
        file,
        sep="\t"
    )

    # ----------------------------------------
    # Clean tool names
    # ----------------------------------------

    df["Tool"] = (

        df["Tool"]

        .str.replace("_sorted","",regex=False)

        .str.replace("_mouse","",regex=False)

        .str.replace("_celegans","",regex=False)

        .str.replace("_human","",regex=False)

        .str.replace(".bed","",regex=False)

    )

    x = range(len(df))

    width = 0.38

    b1 = ax.bar(

        [i-width/2 for i in x],

        df["TPR"],

        width,

        color="#4C72B0",

        label="TPR"

    )

    b2 = ax.bar(

        [i+width/2 for i in x],

        df["FPR"],

        width,

        color="#DD8452",

        label="FPR"

    )

    if legend_handles is None:

        legend_handles = [b1, b2]

        legend_labels = ["TPR", "FPR"]

    ax.set_title(
        title,
        fontsize=13
    )

    ax.set_ylim(
        0,
        1
    )

    ax.set_ylabel(
        "Rate"
    )

    ax.set_xticks(list(x))

    ax.set_xticklabels(

        df["Tool"],

        rotation=45,

        ha="right"

    )

    ax.grid(

        axis="y",

        alpha=0.3

    )

# ======================================================
# Common Legend
# ======================================================

fig.legend(

    legend_handles,

    legend_labels,

    loc="lower center",

    ncol=2,

    frameon=False,

    fontsize=11

)

plt.tight_layout(

    rect=[0,0.05,1,1]

)

# ======================================================
# Save
# ======================================================

plt.savefig(

    OUTPUT_DIR/"Fig7.png",

    dpi=600,

    bbox_inches="tight"

)

plt.savefig(

    OUTPUT_DIR/"Fig7.pdf",

    bbox_inches="tight"

)

plt.close()

print("Done.")