# figure1D_amino_acid_composition.py

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

AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")

# =====================================================
# CHECK INPUT DIRECTORY
# =====================================================

if not INPUT_DIR.exists():
    raise FileNotFoundError(
        f"Input directory does not exist:\n{INPUT_DIR}"
    )

files = sorted(INPUT_DIR.glob("*_aa_count*.tsv"))

if not files:
    raise FileNotFoundError(
        f"No amino acid composition files found in:\n{INPUT_DIR}"
    )

# =====================================================
# LOAD FILES
# =====================================================

all_tools = []

required_columns = {
    "Character",
    "Proportion"
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

    tool = file.stem.replace("_aa_count", "")

    df = df.copy()
    df["Tool"] = tool

    all_tools.append(df)

if not all_tools:
    raise ValueError(
        "No valid amino acid composition files were loaded."
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
    columns="Character",
    values="Proportion",
    aggfunc="sum",
    fill_value=0
)

# Ensure all amino acids exist

plot_df = plot_df.reindex(
    columns=AA_ORDER,
    fill_value=0
)

# Alphabetical ordering of tools

plot_df = plot_df.sort_index()

# =====================================================
# PLOT
# =====================================================

ax = plot_df.plot(
    kind="bar",
    stacked=True,
    figsize=(14, 8)
)

plt.title("Fig 1D")

plt.xlabel("LCR Detection Tool")
plt.ylabel("Proportion")

plt.legend(
    title="Amino Acid",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    ncol=1
)

plt.tight_layout()

output_file = (
    OUTPUT_DIR /
    "Figure1D_AminoAcidComposition.png"
)

plt.savefig(
    output_file,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print(f"\nSaved: {output_file}")