import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# ======================================================
# Input / Output
# ======================================================

INPUT_DIR = Path(sys.argv[1])
OUTPUT_DIR = Path(sys.argv[2])

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

# ======================================================
# Files
# ======================================================

panels = [

    ("gene_length_summary.tsv", "A. Gene Length"),

    ("lcr_count_summary.tsv", "B. LCR Count"),

    ("coverage_summary.tsv", "C. Coverage"),

    ("entropy_ratio_summary.tsv", "D. Entropy Ratio")

]

# ======================================================
# Plot
# ======================================================

fig, axes = plt.subplots(
    2,
    2,
    figsize=(16, 12)
)

axes = axes.flatten()

legend_handles = None
legend_labels = None

for ax, (filename, title) in zip(axes, panels):

    file = INPUT_DIR / filename

    if not file.exists():

        print(f"Skipping {filename}")

        ax.axis("off")

        continue

    df = pd.read_csv(
        file,
        sep="\t"
    )

    required = {
        "Tool",
        "TPR",
        "FPR"
    }

    if not required.issubset(df.columns):

        print(
            f"Skipping {filename}: Missing columns "
            f"{required - set(df.columns)}"
        )

        ax.axis("off")

        continue

    # -------------------------------------------------
    # Clean tool names
    # -------------------------------------------------

    df = df.copy()

    df["Tool"] = (
        df["Tool"]
        .astype(str)
        .str.replace(".bed", "", regex=False)
        .str.replace("_sorted", "", regex=False)
    )

    # -------------------------------------------------
    # Plot
    # -------------------------------------------------

    x = range(len(df))

    width = 0.38

    b1 = ax.bar(
        [i - width / 2 for i in x],
        df["TPR"],
        width,
        color="#4C72B0",
        label="TPR"
    )

    b2 = ax.bar(
        [i + width / 2 for i in x],
        df["FPR"],
        width,
        color="#DD8452",
        label="FPR"
    )

    if legend_handles is None:

        legend_handles = [b1, b2]

        legend_labels = ["TPR", "FPR"]

    ax.set_title(title)

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
    rect=[0, 0.05, 1, 1]
)

# ======================================================
# Save
# ======================================================

output_file = OUTPUT_DIR / "Fig7_summary.png"

plt.savefig(
    output_file,
    dpi=600,
    bbox_inches="tight"
)

plt.close(fig)

print(f"Saved: {output_file}")