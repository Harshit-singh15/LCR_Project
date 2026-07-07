# figure1E_entropy_boxplot.py

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

# =====================================================
# CHECK INPUT DIRECTORY
# =====================================================

if not INPUT_DIR.exists():
    raise FileNotFoundError(
        f"Input directory does not exist:\n{INPUT_DIR}"
    )

files = sorted(INPUT_DIR.glob("*_SNS*.tsv"))

if not files:
    raise FileNotFoundError(
        f"No Shannon entropy files found in:\n{INPUT_DIR}"
    )

# =====================================================
# LOAD FILES
# =====================================================

all_tools = []

required_column = "Shannon_Entropy"

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

    if required_column not in df.columns:
        print(
            f"Skipping {file.name}: "
            f"'{required_column}' column not found."
        )
        continue

    tool = file.stem.replace("_SNS", "")

    df = df.copy()
    df["Tool"] = tool

    all_tools.append(df)

if not all_tools:
    raise ValueError(
        "No valid Shannon entropy files were loaded."
    )

# =====================================================
# COMBINE
# =====================================================

combined = pd.concat(
    all_tools,
    ignore_index=True
)

# Alphabetical ordering

tool_order = sorted(
    combined["Tool"].unique()
)

combined["Tool"] = pd.Categorical(
    combined["Tool"],
    categories=tool_order,
    ordered=True
)

# =====================================================
# PLOT
# =====================================================

plt.figure(figsize=(14, 8))

sns.boxplot(
    data=combined,
    x="Tool",
    y="Shannon_Entropy",
    hue="Tool",
    order=tool_order,
    palette="tab20",
    linewidth=1.2,
    medianprops={
        "color": "black",
        "linewidth": 2
    },
    showfliers=False,
    dodge=False,
    legend=False
)

plt.title("Fig 1E")

plt.xlabel("LCR Detection Tool")
plt.ylabel("Shannon Entropy")

plt.xticks(
    rotation=45,
    ha="right"
)

plt.tight_layout()

output_file = (
    OUTPUT_DIR /
    "Figure1E_Entropy_Boxplot.png"
)

plt.savefig(
    output_file,
    dpi=600,
    bbox_inches="tight"
)

plt.close()

print(f"\nSaved: {output_file}")