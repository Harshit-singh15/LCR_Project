import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# =====================================================
# DIRECTORIES
# =====================================================

INPUT_DIR = Path(
    "Fig2/04_metrics/peptide_counts"
)

OUTPUT_DIR = Path(
    "Fig2/05_plots"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

# =====================================================
# TOP MOTIFS
# =====================================================

TOP_MOTIFS = [

    "L","S","P","E","A",
    "G","K","R","V","T",
    "Q","D","I","C","F",
    "H","N","SS","LL","Y"

]

# =====================================================
# BUILD DATAFRAME
# =====================================================

rows = []

for consensus in range(1,14):

    file = (
        INPUT_DIR /
        f"consensus_{consensus}_peptide_counts.tsv"
    )

    df = pd.read_csv(
        file,
        sep="\t"
    )

    row_dict = {
        motif:0
        for motif in TOP_MOTIFS
    }

    other = 0

    for _, r in df.iterrows():

        motif = r["Best-Peptide"]
        prop = r["Proportion"]

        if motif in TOP_MOTIFS:

            row_dict[motif] += prop

        else:

            other += prop

    row_dict["Other"] = other
    row_dict["Consensus"] = consensus

    rows.append(row_dict)

plot_df = pd.DataFrame(rows)

plot_df = (
    plot_df
    .fillna(0)
    .set_index("Consensus")
)

# =====================================================
# COLORS
# =====================================================

colors = plt.cm.tab20.colors[:20]

plot_colors = list(colors)
plot_colors.append("lightgray")

# =====================================================
# PLOT
# =====================================================

fig, ax = plt.subplots(
    figsize=(14,8)
)

plot_df[
    TOP_MOTIFS + ["Other"]
].plot(

    kind="bar",

    stacked=True,

    ax=ax,

    color=plot_colors,

    width=0.85
)

ax.set_xlabel(
    "Consensus level"
)

ax.set_ylabel(
    "Proportion of LCRs"
)

ax.set_title(
    "Motif composition across consensus levels"
)

ax.legend(

    title="Motif",

    bbox_to_anchor=(1.02,1),

    loc="upper left",

    fontsize=8
)

plt.tight_layout()

plt.savefig(
    OUTPUT_DIR /
    "Fig2A_motif_composition.png",
    dpi=600,
    bbox_inches="tight"
)

plt.show()