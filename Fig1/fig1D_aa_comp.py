# figure1D_amino_acid_composition.py

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

# =====================================================
# CONFIG
# =====================================================

INPUT_DIR = "celegans/outputs_prerequisite/Amino_acid"

OUTPUT_DIR = "celegans/Fig_outputs"

Path(OUTPUT_DIR).mkdir(
    parents=True,
    exist_ok=True
)

AA_ORDER = list(
    "ACDEFGHIKLMNPQRSTVWY"
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
    Path(INPUT_DIR).glob("*_aa_count*")
)

if len(files) == 0:

    raise FileNotFoundError(
        f"No amino acid files found in {INPUT_DIR}"
    )

for file in files:

    print(f"Loading {file.name}")

    tool = file.name.split(
        "_aa_count"
    )[0]

    df = pd.read_csv(
        file,
        sep="\t"
    )

    required = {
        "Character",
        "Proportion"
    }

    if not required.issubset(df.columns):

        print(
            f"Skipping {file.name}"
        )
        continue

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
    columns="Character",
    values="Proportion"
)

plot_df = plot_df.reindex(
    columns=AA_ORDER,
    fill_value=0
)

plot_df = plot_df.reindex(
    tool_order
)

plot_df = plot_df.dropna(
    how="all"
)

plot_df = plot_df.fillna(0)

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

plt.title(
    "C. elegans Figure 1D: Amino Acid Composition"
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