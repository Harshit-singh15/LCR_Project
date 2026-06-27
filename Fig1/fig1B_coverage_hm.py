# figure1B_coverage_heatmap.py

from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# =====================================================
# CONFIG
# =====================================================

INPUT_DIR = "celegans\\outputs_prerequisite\\LCR_Coverage"

OUTPUT_DIR = "celegans\\Fig_outputs"

Path(OUTPUT_DIR).mkdir(
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
    "alcor_mode1_masked_celegans":"AlcoR M1",
    "alcor_mode2_masked_celegans":"AlcoR M2",
    "dotplot_celegans":"Dotplot",
    "flps_default_celegans":"fLPS",
    "flps_strict_celegans":"fLPS Strict",
    "flps2_default_celegans":"fLPS 2.0",
    "flps2_strict_celegans":"fLPS 2.0 Strict",
    "lcrfinder_celegans":"LCRFinder",
    "seg_celegans":"SEG",
    "seg_intermediate_celegans":"SEG Intermediate",
    "seg_strict_celegans":"SEG Strict",
    "treks_combined_celegans":"T-REKS",
    "xstream_m1_celegans":"XSTREAM"
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
            f"Skipping {file.name}"
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

heatmap_df = combined.pivot(
    index="Tool",
    columns="Category",
    values="Count"
)

heatmap_df = heatmap_df.reindex(
    columns=coverage_order,
    fill_value=0
)

heatmap_df = heatmap_df.reindex(
    tool_order
)

heatmap_df = heatmap_df.dropna(
    how="all"
)

heatmap_df = heatmap_df.fillna(0)

heatmap_df.index = [
    tool_labels.get(x, x)
    for x in heatmap_df.index
]

# =====================================================
# PLOT
# =====================================================

plt.figure(
    figsize=(12,8)
)

sns.heatmap(
    heatmap_df,
    cmap="YlOrRd",
    annot=True,
    fmt=".0f",
    linewidths=0.5
)

plt.title(
    "C. elegans Figure 1B: LCR Coverage Distribution"
)

plt.xlabel(
    "Coverage Category (%)"
)

plt.ylabel(
    "LCR Detection Tool"
)

plt.tight_layout()

output_file = (
    Path(OUTPUT_DIR)
    /
    "Figure1B_Coverage_Heatmap.png"
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