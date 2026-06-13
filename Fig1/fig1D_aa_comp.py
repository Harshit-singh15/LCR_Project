# figure1D_amino_acid_composition.py

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

INPUT_DIR = "outputs_prerequisite/Amino_acid"

all_tools = []

for file in Path(INPUT_DIR).glob("*_aa_count*"):

    tool = file.name.split("_aa_count")[0]

    df = pd.read_csv(file, sep=r"\s+")

    df["Tool"] = tool

    all_tools.append(df)

combined = pd.concat(all_tools, ignore_index=True)

plot_df = combined.pivot(
    index="Tool",
    columns="Character",
    values="Proportion"
)

ax = plot_df.plot(
    kind="bar",
    stacked=True,
    figsize=(14,8)
)

plt.title("Figure 1D: Amino Acid Composition")
plt.xlabel("Tool")
plt.ylabel("Proportion")

plt.legend(
    title="Amino Acid",
    bbox_to_anchor=(1.02,1),
    loc="upper left",
    ncol=1
)

plt.tight_layout()

plt.savefig("Fig_outputs/Figure1D_AminoAcidComposition.png", dpi=300)
plt.savefig("Fig_outputs/Figure1D_AminoAcidComposition.pdf")

plt.show()