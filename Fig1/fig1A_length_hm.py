# figure1A_length_heatmap.py

from pathlib import Path
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# =====================================================
# INPUTS
# =====================================================

INPUT_DIR = sys.argv[1]
OUTPUT_DIR = sys.argv[2]

Path(OUTPUT_DIR).mkdir(
    parents=True,
    exist_ok=True
)

length_order = [
    "0-10",
    "10-20",
    "20-50",
    "50-100",
    "100-200",
    "200+"
]

# =====================================================
# LOAD FILES
# =====================================================

all_tools = []

files = sorted(Path(INPUT_DIR).glob("*.tsv"))

if not files:
    raise FileNotFoundError(
        f"No TSV files found in {INPUT_DIR}"
    )

for file in files:

    print(f"Loading {file.name}")

    df = pd.read_csv(
        file,
        sep="\t"
    )

    required = {"Category", "Count"}

    if not required.issubset(df.columns):
        print(
            f"Skipping {file.name} "
            f"(missing required columns)"
        )
        continue

    # Tool name = filename without extension
    tool = file.stem.replace("_categorized", "")

    df["Tool"] = tool

    all_tools.append(df)

if not all_tools:
    raise ValueError("No valid TSV files found.")

# =====================================================
# COMBINE
# =====================================================

combined = pd.concat(
    all_tools,
    ignore_index=True
)

heatmap_df = combined.pivot(
    index="Tool",
    columns="Category",
    values="Count"
)

# Keep category order
heatmap_df = heatmap_df.reindex(
    columns=length_order,
    fill_value=0
)

# Fill missing values
heatmap_df = heatmap_df.fillna(0)

# Sort tools alphabetically
heatmap_df = heatmap_df.sort_index()

# =====================================================
# PLOT
# =====================================================

plt.figure(figsize=(12, 8))

sns.heatmap(
    heatmap_df,
    cmap="YlOrRd",
    annot=True,
    fmt=".0f",
    linewidths=0.5
)

plt.title("Fig 1A")

plt.xlabel("Length Category (aa)")
plt.ylabel("LCR Detection Tool")

plt.tight_layout()

output_file = (
    Path(OUTPUT_DIR)
    / "Figure1A_Length_Heatmap.png"
)

plt.savefig(
    output_file,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print(f"\nSaved: {output_file}")