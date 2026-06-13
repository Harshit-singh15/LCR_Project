# figure1B_coverage_heatmap.py

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

INPUT_DIR = "outputs_prerequisite/LCR_Coverage"

coverage_order = [
    "0-20",
    "20-40",
    "40-60",
    "60-80",
    "80-100"
]

all_tools = []

for file in Path(INPUT_DIR).glob("*.tsv"):

    tool = file.stem.replace("_categorized", "")

    df = pd.read_csv(file, sep="\t")
    df["Tool"] = tool

    all_tools.append(df)

combined = pd.concat(all_tools, ignore_index=True)

heatmap_df = combined.pivot(
    index="Tool",
    columns="Category",
    values="Count"
)

heatmap_df = heatmap_df[coverage_order]

plt.figure(figsize=(10,8))

sns.heatmap(
    heatmap_df,
    cmap="YlOrRd",
    annot=True,
    fmt=".0f"
)

plt.title("Figure 1B: LCR Coverage Distribution")
plt.xlabel("Coverage Category (%)")
plt.ylabel("Tool")

plt.tight_layout()

plt.savefig("Fig_outputs/Figure1B_Coverage_Heatmap.png", dpi=300)
plt.savefig("Fig_outputs/Figure1B_Coverage_Heatmap.pdf")

plt.show() 