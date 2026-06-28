import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# ==========================================================
# INPUT / OUTPUT
# ==========================================================

input_file = r"fruitfly\dataforFig8\fig8B_missing_residue_complexity.tsv"

output_dir = Path(r"fruitfly\Fig_outputs\Fig8")
output_dir.mkdir(
    parents=True,
    exist_ok=True
)

# ==========================================================

print("Loading data...")

df = pd.read_csv(
    input_file,
    sep="\t"
)

print(f"Regions : {len(df):,}")

# ==========================================================
# BINNING
# ==========================================================

bins = np.arange(0, 110, 10)

labels = list(range(10))

df["Mut_Bin"] = pd.cut(

    df["Mutation_Percent"],

    bins=bins,

    labels=labels,

    include_lowest=True,

    right=False

)

df["Freq_Bin"] = pd.cut(

    df["Most_Frequent_AA_Percent"],

    bins=bins,

    labels=labels,

    include_lowest=True,

    right=False

)

# Handle exact value of 100

df.loc[
    df["Mutation_Percent"] == 100,
    "Mut_Bin"
] = 9

df.loc[
    df["Most_Frequent_AA_Percent"] == 100,
    "Freq_Bin"
] = 9

df["Mut_Bin"] = df["Mut_Bin"].astype(int)

df["Freq_Bin"] = df["Freq_Bin"].astype(int)

# ==========================================================
# COUNT REGIONS
# ==========================================================

heat = (

    df.groupby(

        ["Freq_Bin", "Mut_Bin"]

    )

    .size()

    .reset_index(name="Count")

)

pivot = heat.pivot(

    index="Freq_Bin",

    columns="Mut_Bin",

    values="Count"

).fillna(0)

# Always generate a complete 10 × 10 matrix

all_bins = list(range(10))

pivot = pivot.reindex(

    index=all_bins,

    columns=all_bins,

    fill_value=0

)

# Reverse Y-axis to match paper

pivot = pivot.sort_index(
    ascending=False
)

# ==========================================================
# LOG SCALE
# ==========================================================

plot_data = np.log10(
    pivot + 1
)

# ==========================================================
# PLOT
# ==========================================================

sns.set_style("white")

plt.figure(figsize=(8, 6))

ax = sns.heatmap(

    plot_data,

    cmap="Reds",

    linewidths=1,

    linecolor="white",

    square=True,

    cbar_kws={

        "label": "Number of Regions"

    }

)

# ==========================================================
# COLORBAR
# ==========================================================

cbar = ax.collections[0].colorbar

cbar.set_ticks([

    np.log10(1 + 1),

    np.log10(10 + 1),

    np.log10(100 + 1),

    np.log10(1000 + 1)

])

cbar.set_ticklabels([

    "1",

    "10",

    "100",

    "1000"

])

# ==========================================================
# AXES
# ==========================================================

ax.set_xlabel(

    "Mutation Percentage (%)",

    fontsize=12

)

ax.set_ylabel(

    "Most Frequent Amino Acid Percentage (%)",

    fontsize=12

)

ax.set_title(

    "Fig. 8B  Missing Residues from PDB Structures",

    fontsize=14

)

ax.set_xticklabels(

    [f"{i*10}-{i*10+10}" for i in range(10)],

    rotation=45,

    ha="right",

    fontsize=10

)

ax.set_yticklabels(

    [f"{i*10}-{i*10+10}" for i in range(9, -1, -1)],

    rotation=0,

    fontsize=10

)

plt.tight_layout()

# ==========================================================
# SAVE
# ==========================================================

plt.savefig(

    output_dir / "Fig8B_missing_residue_heatmap.png",

    dpi=600,

    bbox_inches="tight"

)

plt.show()