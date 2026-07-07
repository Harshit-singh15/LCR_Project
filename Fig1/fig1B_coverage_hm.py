# figure1B_coverage_heatmap.py

from pathlib import Path
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# =====================================================
# INPUTS
# =====================================================

INPUT_DIR = Path(sys.argv[1])
OUTPUT_DIR = Path(sys.argv[2])

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

coverage_order = [
    "0-20",
    "20-40",
    "40-60",
    "60-80",
    "80-100"
]

# =====================================================
# CHECK INPUT DIRECTORY
# =====================================================

if not INPUT_DIR.exists():
    raise FileNotFoundError(
        f"Input directory does not exist:\n{INPUT_DIR}"
    )

files = sorted(INPUT_DIR.glob("*.tsv"))

if not files:
    raise FileNotFoundError(
        f"No TSV files found in:\n{INPUT_DIR}"
    )

# =====================================================
# LOAD FILES
# =====================================================

all_tools = []

required_columns = {
    "Category",
    "Count"
}

for file in files:

    print(f"Loading {file.name}")

    try:
        df = pd.read_csv(
            file,
            sep="\t"
        )
    except Exception as e:
        print(f"Skipping {file.name}: {e}")
        continue

    if not required_columns.issubset(df.columns):
        print(
            f"Skipping {file.name}: "
            f"Missing columns "
            f"{required_columns - set(df.columns)}"
        )
        continue

    tool = file.stem.replace(
        "_categorized",
        ""
    )

    df = df.copy()
    df["Tool"] = tool

    all_tools.append(df)

if not all_tools:
    raise ValueError(
        "No valid TSV files were found."
    )

# =====================================================
# COMBINE
# =====================================================

combined = pd.concat(
    all_tools,
    ignore_index=True
)

heatmap_df = combined.pivot_table(
    index="Tool",
    columns="Category",
    values="Count",
    aggfunc="sum",
    fill_value=0
)

# Ensure category order

heatmap_df = heatmap_df.reindex(
    columns=coverage_order,
    fill_value=0
)

# Alphabetical order of tools

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

plt.title("Fig 1B")

plt.xlabel("Coverage Category (%)")
plt.ylabel("LCR Detection Tool")

plt.tight_layout()

output_file = OUTPUT_DIR / "Figure1B_Coverage_Heatmap.png"

plt.savefig(
    output_file,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print(f"\nSaved: {output_file}")