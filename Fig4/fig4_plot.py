import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import sys

# ======================================================
# Input / Output
# ======================================================

INPUT_FILE = Path(sys.argv[1])

OUTPUT_DIR = Path(sys.argv[2])

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

# ======================================================
# Read matrix
# ======================================================

df = pd.read_csv(
    INPUT_FILE,
    sep="\t",
    index_col=0
)

# ======================================================
# Clean method names
# ======================================================

def clean_name(name):

    name = str(name)

    name = name.replace(".bed", "")
    name = name.replace("_sorted", "")

    return name

df.index = [clean_name(x) for x in df.index]
df.columns = [clean_name(x) for x in df.columns]

# ======================================================
# Lower triangle only
# ======================================================

mask = np.triu(
    np.ones_like(df, dtype=bool),
    k=1
)

# ======================================================
# Plot
# ======================================================

plt.figure(
    figsize=(10, 8)
)

sns.heatmap(

    df,

    mask=mask,

    cmap="Reds",

    vmin=0,
    vmax=1,

    annot=True,
    fmt=".2f",

    annot_kws={
        "size": 8
    },

    square=True,

    linewidths=0.5,

    linecolor="white",

    cbar_kws={
        "label": "Jaccard similarity"
    }

)

plt.title(
    "Fig 4"
)

plt.xticks(
    rotation=45,
    ha="right"
)

plt.yticks(
    rotation=0
)

plt.tight_layout()

# ======================================================
# Save
# ======================================================

output_file = OUTPUT_DIR / "Fig4_Jaccard.png"

plt.savefig(
    output_file,
    dpi=600,
    bbox_inches="tight"
)

plt.close()

print(f"Saved: {output_file}")