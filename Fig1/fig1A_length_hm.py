# figure1A_length_heatmap.py

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

INPUT_DIR = "outputs_prerequisite/LCR_Length"

length_order = [
    "0-10",
    "10-20",
    "20-50",
    "50-100",
    "100-200",
    "200+"
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

heatmap_df = heatmap_df[length_order]

plt.figure(figsize=(10, 8))

sns.heatmap(
    heatmap_df,
    cmap="YlOrRd",
    annot=True,
    fmt=".0f"
)

plt.title("Figure 1A: LCR Length Distribution")
plt.xlabel("Length Category (aa)")
plt.ylabel("Tool")

plt.tight_layout()

plt.savefig("Fig_outputs/Figure1A_Length_Heatmap.png", dpi=300)
plt.savefig("Fig_outputs/Figure1A_Length_Heatmap.pdf")

plt.show()