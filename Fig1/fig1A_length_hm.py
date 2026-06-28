# figure1A_length_heatmap.py

from pathlib import Path

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# =====================================================
# CONFIG
# =====================================================

INPUT_DIR = r"arabidopsis\dataforFig1\LCR_Length"

OUTPUT_DIR = r"arabidopsis\Fig_outputs\Fig1"

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

heatmap_df = combined.pivot(
    index="Tool",
    columns="Category",
    values="Count"
)

# enforce category order

heatmap_df = heatmap_df.reindex(
    columns=length_order,
    fill_value=0
)

# enforce tool order

heatmap_df = heatmap_df.reindex(
    tool_order
)
heatmap_df.index = [
    tool_labels.get(x, x)
    for x in heatmap_df.index
]

# remove tools not present

heatmap_df = heatmap_df.dropna(
    how="all"
)

# replace remaining NaN values

heatmap_df = heatmap_df.fillna(0)

# =====================================================
# PLOT
# =====================================================

plt.figure(
    figsize=(12, 8)
)

sns.heatmap(
    heatmap_df,
    cmap="YlOrRd",
    annot=True,
    fmt=".0f",
    linewidths=0.5
)

organism = files[0].stem.split("_")[-2]

plt.title(
    f"{organism.capitalize()} Fig 1A: LCR Length Distribution"
)

plt.xlabel(
    "Length Category (aa)"
)

plt.ylabel(
    "LCR Detection Tool"
)

plt.tight_layout()

output_file = (
    Path(OUTPUT_DIR)
    /
    "Figure1A_Length_Heatmap.png"
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
