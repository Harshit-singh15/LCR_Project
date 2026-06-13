import pandas as pd
import matplotlib.pyplot as plt

files = [
    ("gene_length_summary.tsv",
     "A. Gene Length"),

    ("lcr_count_summary.tsv",
     "B. LCR Count"),

    ("coverage_summary.tsv",
     "C. Coverage"),

    ("entropy_ratio_summary.tsv",
     "D. Entropy Ratio")
]

fig, axes = plt.subplots(
    2,
    2,
    figsize=(18,12)
)

axes = axes.flatten()

for ax, (file, title) in zip(axes, files):

    df = pd.read_csv(
        f"fig7/{file}",
        sep="\t"
    )

    x = range(len(df))

    width = 0.4

    ax.bar(
        [i-width/2 for i in x],
        df["TPR"],
        width,
        label="TPR"
    )

    ax.bar(
        [i+width/2 for i in x],
        df["FPR"],
        width,
        label="FPR"
    )

    ax.set_xticks(x)

    ax.set_xticklabels(
        df["Tool"],
        rotation=90
    )

    ax.set_ylim(0,1)

    ax.set_title(title)

    ax.legend()

plt.tight_layout()

plt.savefig(
    "Fig7/Figure7.png",
    dpi=300
)

plt.show()