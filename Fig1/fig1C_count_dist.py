# figure1C_count_distribution.py

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

# =====================================================
# CONFIG
# =====================================================

INPUT_DIR = "Celegans/outputs_prerequisite/LCR_Count"

OUTPUT_DIR = "Celegans/Fig_outputs"

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
    "alcor_mode1_masked_celegans",
    "alcor_mode2_masked_celegans",
    "dotplot_celegans",
    "flps_default_celegans",
    "flps_strict_celegans",
    "flps2_default_celegans",
    "flps2_strict_celegans",
    "lcrfinder_celegans",
    "seg_celegans",
    "seg_intermediate_celegans",
    "seg_strict_celegans",
    "treks_combined_celegans",
    "xstream_m1_celegans"
]

tool_labels = {
    "alcor_mode1_masked_celegans": "AlcoR M1",
    "alcor_mode2_masked_celegans": "AlcoR M2",
    "dotplot_celegans": "Dotplot",
    "flps_default_celegans": "fLPS",
    "flps_strict_celegans": "fLPS Strict",
    "flps2_default_celegans": "fLPS 2.0",
    "flps2_strict_celegans": "fLPS 2.0 Strict",
    "lcrfinder_celegans": "LCRFinder",
    "seg_celegans": "SEG",
    "seg_intermediate_celegans": "SEG Intermediate",
    "seg_strict_celegans": "SEG Strict",
    "treks_combined_celegans": "T-REKS",
    "xstream_m1_celegans": "XSTREAM"
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

    tool = file.stem.replace(
        "_categorized",
        ""
    )

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

plt.title(
    "C. elegans Figure 1C: Number of LCRs per Protein"
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