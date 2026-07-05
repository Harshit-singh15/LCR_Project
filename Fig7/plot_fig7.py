import re
from pathlib import Path
import sys
import pandas as pd
import matplotlib.pyplot as plt

# ======================================================
# Input / Output
# ======================================================

INPUT_DIR = Path(
    sys.argv[1]  # Input directory containing summary tables
)

OUTPUT_DIR = Path(
    sys.argv[2]  # Output directory for combined figure
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
    # Clean tool names (organism-independent)
    # ----------------------------------------

    tool_labels = {
        "alcor_mode1_masked": "AlcoR M1",
        "alcor_mode2_masked": "AlcoR M2",
        "dotplot": "Dotplot",
        "flps_default": "fLPS",
        "flps_strict": "fLPS Strict",
        "flps2_default": "fLPS 2.0",
        "flps2_strict": "fLPS 2.0 Strict",
        "lcrfinder": "LCRFinder",
        "seg": "SEG",
        "seg_intermediate": "SEG Intermediate",
        "seg_strict": "SEG Strict",
        "treks_combined": "T-REKS",
        "xstream_m1": "XSTREAM"
    }

    organisms = {
        "human",
        "mouse",
        "zebrafish",
        "celegans",
        "Fruitfly",
        "arabidopsis",
        "ecoli",
        "yeast"
    }

    clean_names = []

    for tool in df["Tool"]:

        tool = tool.replace(".bed", "")
        tool = tool.replace("_sorted", "")

        for org in organisms:
            suffix = "_" + org
            if tool.endswith(suffix):
                tool = tool[:-len(suffix)]
                break

        clean_names.append(
            tool_labels.get(tool, tool)
        )

    df["Tool"] = clean_names

    # ----------------------------------------
    # Plot bars
    # ----------------------------------------

    x = range(len(df))
    width = 0.38

    b1 = ax.bar(
        [i - width/2 for i in x],
        df["TPR"],
        width,
        color="#4C72B0",
        label="TPR"
    )

    b2 = ax.bar(
        [i + width/2 for i in x],
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

    ax.set_ylim(0, 1)

    ax.set_ylabel("Rate")

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

    OUTPUT_DIR/"Fig7_summary.png",

    dpi=600,

    bbox_inches="tight"

)


plt.close()

print("Done.")