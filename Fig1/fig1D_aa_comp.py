# figure1D_amino_acid_composition.py

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

# =====================================================
# CONFIG
# =====================================================

INPUT_DIR = r"arabidopsis\dataforFig1\Amino_acid"

OUTPUT_DIR = r"arabidopsis\Fig_outputs\Fig1"

Path(OUTPUT_DIR).mkdir(
    parents=True,
    exist_ok=True
)

AA_ORDER = list(
    "ACDEFGHIKLMNPQRSTVWY"
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
    "alcor_mode1_masked":"AlcoR M1",
    "alcor_mode2_masked":"AlcoR M2",
    "dotplot":"Dotplot",
    "flps_default":"fLPS",
    "flps_strict":"fLPS Strict",
    "flps2_default":"fLPS 2.0",
    "flps2_strict":"fLPS 2.0 Strict",
    "lcrfinder":"LCRFinder",
    "seg":"SEG",
    "seg_intermediate":"SEG Intermediate",
    "seg_strict":"SEG Strict",
    "treks_clustalw":"T-REKS",
    "xstream_m1":"XSTREAM"
}

# =====================================================
# LOAD FILES
# =====================================================

all_tools = []

files = sorted(Path(INPUT_DIR).glob("*_aa_count*"))

if not files:
    raise FileNotFoundError(
        f"No amino acid files found in {INPUT_DIR}"
    )

# Known organisms (extend if required)
organisms = [
    "human",
    "mouse",
    "zebrafish",
    "celegans",
    "arabidopsis",
    "Fruitfly",
    "ecoli",
    "yeast"
]

for file in files:

    print(f"Loading {file.name}")

    # Example:
    # dotplot_zebrafish_aa_count.tsv
    # flps_default_mouse_aa_count.tsv

    tool = file.stem

    # remove suffix
    tool = tool.replace("_aa_count", "")

    # remove organism if present
    for org in organisms:
        suffix = "_" + org
        if tool.endswith(suffix):
            tool = tool[:-len(suffix)]
            break

    print("Detected tool:", tool)

    df = pd.read_csv(
        file,
        sep="\t"
    )

    required = {"Character", "Proportion"}

    if not required.issubset(df.columns):
        print(
            f"Skipping {file.name} "
            f"(missing columns: {required - set(df.columns)})"
        )
        continue

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

plot_df = combined.pivot(
    index="Tool",
    columns="Character",
    values="Proportion"
)

plot_df = plot_df.reindex(
    columns=AA_ORDER,
    fill_value=0
)

# Keep only tools that actually exist
available_tools = [
    t for t in tool_order
    if t in plot_df.index
]

missing_tools = sorted(
    set(tool_order) - set(available_tools)
)

if missing_tools:
    print("\nMissing tools:")
    print(missing_tools)

plot_df = plot_df.reindex(available_tools)
plot_df = plot_df.fillna(0)

if plot_df.empty:
    raise ValueError(
        "No matching tools found after parsing filenames."
    )

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
    figsize=(14, 8)
)

organism = Path(INPUT_DIR).parent.parent.name

plt.title(
    f"{organism.capitalize()} Fig 1D: Amino Acid Composition"
)

plt.xlabel(
    "LCR Detection Tool"
)

plt.ylabel(
    "Proportion"
)

plt.legend(
    title="Amino Acid",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    ncol=1
)

plt.tight_layout()

plt.savefig(
    Path(OUTPUT_DIR) /
    "Figure1D_AminoAcidComposition.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()