import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from pathlib import Path
import re

# ============================================
# CONFIGURATION
# ============================================
COMPLEXITY_DIR = r"SupplyFig\complexity"
BIN_SIZE = 2

# --------------------------------------------
# Find and Validate Complexity Files
# --------------------------------------------
files = sorted(Path(COMPLEXITY_DIR).glob("*complexity*"))
if len(files) == 0:
    raise FileNotFoundError("No complexity files found.")

def get_k(filepath):
    name = filepath.stem
    m = re.search(r"(\d+)", name)
    if not m:
        raise ValueError(f"Cannot determine k from {name}")
    return int(m.group(1))

# --------------------------------------------
# Read and Concat Data
# --------------------------------------------
dfs = []
for f in files:
    k = get_k(f)
    df_temp = pd.read_csv(f, sep="\t")
    df_temp["k"] = k
    dfs.append(df_temp)

df = pd.concat(dfs, ignore_index=True)

# Clean numeric columns
df["Mutation_Percent"] = pd.to_numeric(df["Mutation_Percent"], errors="coerce")
df["Most_Frequent_AA_Percent"] = pd.to_numeric(df["Most_Frequent_AA_Percent"], errors="coerce")
df = df.dropna(subset=["Mutation_Percent", "Most_Frequent_AA_Percent", "k"])

# --------------------------------------------
# Bin Data (2D grid allocation)
# --------------------------------------------
df["x_bin"] = np.floor(df["Mutation_Percent"] / BIN_SIZE) * BIN_SIZE
df["y_bin"] = np.floor(df["Most_Frequent_AA_Percent"] / BIN_SIZE) * BIN_SIZE

binned = df.groupby(["k", "x_bin", "y_bin"]).size().reset_index(name="n")
global_max = binned["n"].max()

# --------------------------------------------
# Structured Grid Plotting (Clean Style)
# --------------------------------------------
fig, axes = plt.subplots(
    4, 4, 
    figsize=(14, 13), 
    sharex=True, 
    sharey=True
)

norm = PowerNorm(gamma=0.5, vmin=0, vmax=global_max)
im = None

# Loop over all 16 grid positions
for idx, ax in enumerate(axes.flatten()):
    k = idx + 1  # 1-indexed for k values
    
    # Clean up axes spines to match target style (Despining)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(False)  # Ensure no background grid lines
    ax.set_facecolor('white') # Clean white background
    
    if k <= 13:
        subset = binned[binned["k"] == k]
        heat = np.zeros((51, 51))

        for _, row in subset.iterrows():
            x = int(row["x_bin"] / BIN_SIZE)
            y = int(row["y_bin"] / BIN_SIZE)
            if 0 <= x < 51 and 0 <= y < 51:
                heat[y, x] = row["n"]

        # Mask 0 values so they show up as pure background white rather than the lowest cmap color
        heat_masked = np.ma.masked_where(heat == 0, heat)

        im = ax.imshow(
            heat_masked,
            origin="lower",
            extent=[0, 100, 0, 100],
            cmap="plasma",
            norm=norm,
            interpolation="nearest",
            aspect="equal"
        )
        
        # Clean title matching your image style (just the raw k-integer)
        ax.set_title(str(k), fontsize=14, pad=8)
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 100)
        
        # Subtle tick mark adjustments
        ax.tick_params(axis='both', which='major', labelsize=11)
        
    else:
        # Completely hide the unused 14, 15, and 16 panels
        ax.axis("off")

# --------------------------------------------
# Labels and Titles
# --------------------------------------------
fig.supxlabel("Mutation (%)", fontsize=16, y=0.04)
fig.supylabel("Most frequent amino acid (%)", fontsize=16, x=0.04)

# Dedicated spatial formatting adjustments to prevent colorbar squeezing
fig.subplots_adjust(right=0.85, left=0.08, bottom=0.1, top=0.92, wspace=0.2, hspace=0.35)

# --------------------------------------------
# Clean Minimalist Colorbar
# --------------------------------------------
colorbar_axis = fig.add_axes([0.89, 0.35, 0.02, 0.35]) # Sized and centered on the right side
cbar = fig.colorbar(im, cax=colorbar_axis)
cbar.set_label("Count", fontsize=14, labelpad=10, rotation=270, va='bottom')
cbar.ax.tick_params(labelsize=11)
cbar.outline.set_visible(False) # Removes the rigid bounding box around colorbar

# Save the styled figure
plt.savefig(
    r"SupplyFig\plots\Suppl_Fig4_LC_small_multiples.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()