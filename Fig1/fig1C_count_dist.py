# figure1C_count_distribution.py

from pathlib import Path
import sys
import pandas as pd
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

count_order = [
    "0",
    "1-5",
    "6-10",
    "11-15",
    "16+"
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
            f"Skipping {file.name}: Missing columns "
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

plot_df = combined.pivot_table(
    index="Tool",
    columns="Category",
    values="Count",
    aggfunc="sum",
    fill_value=0
)

# Keep category order

plot_df = plot_df.reindex(
    columns=count_order,
    fill_value=0
)

# Alphabetical order of tools

plot_df = plot_df.sort_index()

# =====================================================
# PLOT
# =====================================================

ax = plot_df.plot(
    kind="bar",
    stacked=True,
    figsize=(12, 8)
)

plt.title("Fig 1C")

plt.xlabel("LCR Detection Tool")
plt.ylabel("Protein Count")

plt.legend(
    title="LCR Count Category",
    bbox_to_anchor=(1.02, 1),
    loc="upper left"
)

plt.tight_layout()

output_file = OUTPUT_DIR / "Figure1C_CountDistribution.png"

plt.savefig(
    output_file,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print(f"\nSaved: {output_file}")