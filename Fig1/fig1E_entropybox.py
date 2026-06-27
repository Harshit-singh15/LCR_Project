# figure1E_entropy_boxplot.py

from pathlib import Path

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# =====================================================
# CONFIG
# =====================================================

INPUT_DIR = "Celegans\\outputs_prerequisite\\ShanonEntropy"

OUTPUT_DIR = "Celegans\\Fig_outputs"

Path(OUTPUT_DIR).mkdir(
    parents=True,
    exist_ok=True
)

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
    Path(INPUT_DIR).glob("*_SNS*")
)

if len(files) == 0:

    raise FileNotFoundError(
        f"No entropy files found in {INPUT_DIR}"
    )

for file in files:

    print(f"Loading {file.name}")

    tool = file.name.replace(
        "_SNS.tsv",
        ""
    ).replace(
        "_SNS",
        ""
    )

    try:

        df = pd.read_csv(
            file,
            sep="\t"
        )

    except Exception as e:

        print(
            f"Skipping {file.name}: {e}"
        )

        continue

    if "Shannon_Entropy" not in df.columns:

        print(
            f"Skipping {file.name}: "
            f"Shannon_Entropy column missing"
        )

        continue

    df["Tool"] = tool

    all_tools.append(df)

if len(all_tools) == 0:

    raise ValueError(
        "No valid entropy files loaded."
    )

# =====================================================
# COMBINE
# =====================================================

combined = pd.concat(
    all_tools,
    ignore_index=True
)

# keep only tools present

tool_order_present = [
    tool
    for tool in tool_order
    if tool in combined["Tool"].unique()
]

combined["Tool"] = pd.Categorical(
    combined["Tool"],
    categories=tool_order_present,
    ordered=True
)

combined["Tool_Label"] = combined[
    "Tool"
].map(tool_labels)

# =====================================================
# PLOT
# =====================================================

plt.figure(
    figsize=(14, 8)
)

sns.boxplot(
    data=combined,
    x="Tool_Label",
    y="Shannon_Entropy",
    hue="Tool_Label",
    order=[
        tool_labels[t]
        for t in tool_order_present
    ],
    palette="tab20",
    linewidth=1.2,
    medianprops={
        "color": "black",
        "linewidth": 2
    },
    showfliers=False,
    legend=False
)

plt.title(
    "C. elegans Figure 1E: Shannon Entropy Distribution"
)

plt.xlabel(
    "LCR Detection Tool"
)

plt.ylabel(
    "Shannon Entropy"
)

plt.xticks(
    rotation=45,
    ha="right"
)

plt.tight_layout()

plt.savefig(
    Path(OUTPUT_DIR) /
    "Figure1E_Entropy_Boxplot.png",
    dpi=600,
    bbox_inches="tight"
)

plt.savefig(
    Path(OUTPUT_DIR) /
    "Figure1E_Entropy_Boxplot.pdf",
    bbox_inches="tight"
)

print(
    "\nSaved Figure1E_Entropy_Boxplot.png"
)

plt.show()