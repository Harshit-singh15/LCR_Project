# figure1E_entropy_boxplot.py

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

INPUT_DIR = "outputs_prerequisite/Diversity"

all_tools = []

for file in Path(INPUT_DIR).iterdir():

    if not file.is_file():
        continue

    tool = file.name.replace("_SNS", "")

    df = pd.read_csv(file)

    df["Tool"] = tool

    all_tools.append(df)

combined = pd.concat(all_tools, ignore_index=True)

# optional: consistent tool ordering
tool_order = sorted(combined["Tool"].unique())

plt.figure(figsize=(14, 8))

# Updated to prevent the v0.14.0 deprecation warning
sns.boxplot(
    data=combined,
    x="Tool",
    y="Shannon_Entropy",
    hue="Tool",            # Added to explicitly map the palette colors to the tools
    legend=False,          # Added to prevent a redundant legend from rendering
    order=tool_order,
    palette="tab20",
    linewidth=1.2,
    medianprops=dict(color="black", linewidth=2),
    showfliers=False
)

plt.xticks(rotation=45, ha="right")

plt.title("Figure 1E: Shannon Entropy Distribution")
plt.xlabel("Tool")
plt.ylabel("Shannon Entropy")

plt.tight_layout()

plt.savefig("Fig_outputs/Figure1E_Entropy_Boxplot.png", dpi=600)
plt.savefig("Fig_outputs/Figure1E_Entropy_Boxplot.pdf")

plt.show()