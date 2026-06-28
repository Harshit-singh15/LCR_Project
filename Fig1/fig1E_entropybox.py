# figure1E_entropy_boxplot.py

from pathlib import Path

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# =====================================================
# CONFIG
# =====================================================

INPUT_DIR = r"arabidopsis\dataforFig1\ShanonEntropy"

OUTPUT_DIR = r"arabidopsis\Fig_outputs\Fig1"

Path(OUTPUT_DIR).mkdir(
    parents=True,
    exist_ok=True
)

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
    Path(INPUT_DIR).glob("*_SNS*")
)

if len(files) == 0:

    raise FileNotFoundError(
        f"No entropy files found in {INPUT_DIR}"
    )

for file in files:

    print(f"Loading {file.name}")

    # -------------------------------------------------
# Robust tool name extraction
# -------------------------------------------------

    tool = file.stem.replace("_SNS", "")

    organisms = {
    "human",
    "mouse",
    "zebrafish",
    "arabidopsis",
    "celegans",
    "ecoli",
    "Fruitfly",
    "yeast"
    }

    for org in organisms:
        suffix = "_" + org
        if tool.endswith(suffix):
            tool = tool[:-len(suffix)]
            break

    print(f"Detected tool : {tool}")

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

available = set(combined["Tool"].unique())

tool_order_present = [
    t for t in tool_order
    if t in available
]

missing = sorted(
    set(tool_order) - available
)

if missing:
    print("\nMissing tools:")
    print(", ".join(missing))

if len(tool_order_present) == 0:
    raise ValueError(
        "No recognised tools found."
    )

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

organism = Path(INPUT_DIR).parent.parent.name

pretty_name = {
    "celegans": "C. elegans",
    "mouse": "Mouse",
    "human": "Human",
    "zebrafish": "Zebrafish",
    "arabidopsis": "Arabidopsis",
    "Fruitfly": "Fruit Fly"
}.get(
    organism,
    organism.capitalize()
)

plt.title(
    f"{pretty_name} Fig 1E: Shannon Entropy Distribution"
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

print(
    "\nSaved Figure1E_Entropy_Boxplot.png"
)

plt.show()