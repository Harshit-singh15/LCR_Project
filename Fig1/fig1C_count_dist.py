# figure1C_count_distribution.py

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

# =====================================================
# CONFIG
# =====================================================

INPUT_DIR = r"ecoli\dataforFig1\LCR_Count"

OUTPUT_DIR = r"ecoli\Fig_outputs\Fig1"

Path(OUTPUT_DIR).mkdir(
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

tool_order = [
    "alcor_mode1_masked",
    "alcor_mode2_masked",
    "dotplot",
    "flps_default",
    "flps_strict",
    "flps2_default",
    "flps2_strict",
    "lcrfinder",
    "seg",
    "seg_intermediate",
    "seg_strict",
    "treks_clustalw",
    "xstream_m1"
]

tool_labels = {
    "alcor_mode1_masked": "AlcoR M1",
    "alcor_mode2_masked": "AlcoR M2",
    "dotplot": "Dotplot",
    "flps_default": "fLPS",
    "flps_strict": "fLPS Strict",
    "flps2_default": "fLPS 2.0",
    "flps2_strict": "fLPS 2.0 Strict",
    "lcrfinder": "LCRFinder",
    "seg": "SEG",
    "seg_intermediate": "SEG Intermediate",
    "seg_strict": "SEG Strict",
    "treks_clustalw": "T-REKS",
    "xstream_m1": "XSTREAM"
}

# =====================================================
# LOAD FILES
# =====================================================

all_tools = []

files = sorted(
    Path(INPUT_DIR).glob("*.tsv")
)

if len(files) == 0:

    raise FileNotFoundError(
        f"No TSV files found in {INPUT_DIR}"
    )

for file in files:

    print(f"Loading {file.name}")

    df = pd.read_csv(
        file,
        sep="\t"
    )

    required = {
        "Category",
        "Count"
    }

    if not required.issubset(df.columns):

        print(
            f"Skipping {file.name} "
            f"(missing required columns)"
        )
        continue

    tool = file.stem.replace("_categorized", "")

    # remove organism suffix automatically
    parts = tool.split("_")

    for i in range(len(parts), 0, -1):
        candidate = "_".join(parts[:i])
        if candidate in tool_order:
            tool = candidate
            break


    df["Tool"] = tool

    all_tools.append(df)

# =====================================================
# COMBINE
# =====================================================

combined = pd.concat(
    all_tools,
    ignore_index=True
)

plot_df = combined.pivot(
    index="Tool",
    columns="Category",
    values="Count"
)

# Ensure category order

plot_df = plot_df.reindex(
    columns=count_order,
    fill_value=0
)

# Ensure tool order

plot_df = plot_df.reindex(
    tool_order
)

# Remove tools not present

plot_df = plot_df.dropna(
    how="all"
)

# Fill remaining NaN

plot_df = plot_df.fillna(0)

# Replace long names with display names

plot_df.index = [
    tool_labels.get(x, x)
    for x in plot_df.index
]

# =====================================================
# PLOT
# =====================================================

ax = plot_df.plot(
    kind="bar",
    stacked=True,
    figsize=(12, 8)
)

organism = files[0].stem.split("_")[-2]

plt.title(
    f"{organism.capitalize()} Fig 1C: No. of LCRs per Protein"
)


plt.xlabel(
    "LCR Detection Tool"
)

plt.ylabel(
    "Protein Count"
)

plt.legend(
    title="LCR Count Category",
    bbox_to_anchor=(1.02, 1),
    loc="upper left"
)

plt.tight_layout()

output_file = (
    Path(OUTPUT_DIR)
    /
    "Figure1C_CountDistribution.png"
)

plt.savefig(
    output_file,
    dpi=300,
    bbox_inches="tight"
)

print(
    f"\nSaved: {output_file}"
)

plt.show()