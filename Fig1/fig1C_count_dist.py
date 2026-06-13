# figure1C_count_distribution.py

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

INPUT_DIR = "outputs_prerequisite/LCR_Count"

count_order = [
    "0",
    "1-5",
    "6-10",
    "11-15",
    "16+"
]

all_tools = []

for file in Path(INPUT_DIR).glob("*.tsv"):

    tool = file.stem.replace("_categorized", "")

    df = pd.read_csv(file, sep="\t")

    df["Tool"] = tool

    all_tools.append(df)

combined = pd.concat(all_tools, ignore_index=True)

plot_df = combined.pivot(
    index="Tool",
    columns="Category",
    values="Count"
)

plot_df = plot_df[count_order]

ax = plot_df.plot(
    kind="bar",
    stacked=True,
    figsize=(12,8)
)

plt.title("Figure 1C: Number of LCRs per Protein")
plt.xlabel("Tool")
plt.ylabel("Protein Count")

plt.legend(
    title="LCR Count Category",
    bbox_to_anchor=(1.02,1),
    loc="upper left"
)

plt.tight_layout()

plt.savefig("Fig_outputs/Figure1C_CountDistribution.png", dpi=300)
plt.savefig("Fig_outputs/Figure1C_CountDistribution.pdf")

plt.show()