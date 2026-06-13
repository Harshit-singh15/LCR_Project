import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# ======================================
# INPUT
# ======================================

input_file = "Fig8\\fig8B_missing_residue_complexity.tsv"

df = pd.read_csv(
    input_file,
    sep="\t"
)

print("Regions:", len(df))

# ======================================
# BINNING
# ======================================

df["Mut_Bin"] = (
    df["Mutation_Percent"] // 10
).astype(int)

df["Freq_Bin"] = (
    df["Most_Frequent_AA_Percent"] // 10
).astype(int)

# Keep inside 0–9
df["Mut_Bin"] = df["Mut_Bin"].clip(0, 9)
df["Freq_Bin"] = df["Freq_Bin"].clip(0, 9)

# ======================================
# COUNT REGIONS PER BIN
# ======================================

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

# FIX: Force the pivot table to have all bins from 0 to 9 
# This ensures there are always exactly 10 rows and 10 columns.
all_bins = list(range(10))
pivot = pivot.reindex(index=all_bins, columns=all_bins, fill_value=0)

# ======================================
# REVERSE Y-AXIS ORDER
# (matches paper orientation)
# ======================================

pivot = pivot.sort_index(ascending=False)

# ======================================
# LOG SCALE
# ======================================

plot_data = np.log10(pivot + 1)

# ======================================
# PLOT
# ======================================

sns.set_style("white")

plt.figure(figsize=(8, 6))

ax = sns.heatmap(
    plot_data,
    cmap="Reds",
    linewidths=1,
    linecolor="white",
    square=True,
    cbar_kws={
        "label": "Count"
    }
)

# ======================================
# COLORBAR LABELS
# ======================================

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

# ======================================
# AXES
# ======================================

ax.set_xlabel(
    "Mutation Percent (binned)"
)

ax.set_ylabel(
    "Most Frequent AA Percent (binned)"
)

# Show actual bin labels (Now safe because dimensions are guaranteed 10x10)
ax.set_xticklabels(
    [str(i) for i in range(10)],
    rotation=0
)

ax.set_yticklabels(
    [str(i) for i in range(9, -1, -1)],
    rotation=0
)

plt.tight_layout()

plt.savefig(
    "Fig8\\Fig8B_missing_residue_heatmap.png",
    dpi=600,
    bbox_inches="tight"
)

plt.savefig(
    "Fig8\\Fig8B_missing_residue_heatmap.pdf",
    bbox_inches="tight"
)

plt.show()