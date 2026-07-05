import re
from pathlib import Path
import sys  
import pandas as pd
import matplotlib.pyplot as plt

# ======================================================
# Input / Output
# ======================================================

INPUT_DIR = Path(sys.argv[1])     # peptide counts
OUTPUT_DIR = Path(sys.argv[2])    # Fig2

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ======================================================
# Read all consensus files
# ======================================================

files = sorted(
    INPUT_DIR.glob("consensus_*_peptide_counts.tsv"),
    key=lambda x: int(re.search(r"consensus_(\d+)", x.name).group(1))
)

plot_rows = []

for file in files:

    consensus = int(
        re.search(
            r"consensus_(\d+)",
            file.name
        ).group(1)
    )

    df = pd.read_csv(file, sep="\t")

    df = df.sort_values(
        "Count",
        ascending=False
    )

    top20 = df.head(20).copy()

    others = df.iloc[20:]

    if not others.empty:

        top20.loc[len(top20)] = {
            "Best-Peptide": "Others",
            "Best-Type": "Other",
            "Count": others["Count"].sum(),
            "Proportion": others["Proportion"].sum()
        }

    top20["Consensus"] = consensus

    plot_rows.append(
        top20[
            [
                "Consensus",
                "Best-Peptide",
                "Best-Type",
                "Proportion"
            ]
        ]
    )

plot_df = pd.concat(
    plot_rows,
    ignore_index=True
)

# ======================================================
# Store motif type
# ======================================================

motif_type = (
    plot_df
    .drop_duplicates("Best-Peptide")
    .set_index("Best-Peptide")["Best-Type"]
    .to_dict()
)

# ======================================================
# Pivot
# ======================================================

pivot = plot_df.pivot(
    index="Consensus",
    columns="Best-Peptide",
    values="Proportion"
).fillna(0)

# Sort motifs by total abundance

totals = pivot.sum(axis=0).sort_values(
    ascending=False
)

cols = list(totals.index)

if "Others" in cols:
    cols.remove("Others")
    cols.append("Others")

pivot = pivot[cols]

# ======================================================
# Hatch patterns
# ======================================================

hatches = {
    "Mono": "",
    "Di": "///",
    "Tri": "xxx",
    "Other": "..."
}

# ======================================================
# Plot
# ======================================================

fig, ax = plt.subplots(
    figsize=(12,7)
)

bottom = None

for peptide in pivot.columns:

    ptype = motif_type.get(
        peptide,
        "Other"
    )

    hatch = hatches.get(
        ptype,
        ""
    )

    ax.bar(
        pivot.index.astype(str),
        pivot[peptide],
        bottom=bottom,
        label=peptide,
        hatch=hatch,
        edgecolor="black",
        linewidth=0.25
    )

    if bottom is None:
        bottom = pivot[peptide].values.copy()
    else:
        bottom += pivot[peptide].values

# ======================================================
# Formatting
# ======================================================

ax.set_xlabel(
    "Consensus Level",
    fontsize=13
)

ax.set_ylabel(
    "Proportion",
    fontsize=13
)

ax.set_title(
    "E. coli: Top peptide motifs across consensus levels",
    fontsize=15
)
        
ax.set_ylim(0,1)

ax.set_xticks(range(len(pivot.index)))
ax.set_xticklabels(
    pivot.index.astype(str)
)

# ======================================================
# Legend
# ======================================================

ax.legend(
    title="Motif",
    bbox_to_anchor=(1.02,1),
    loc="upper left",
    ncol=2,
    fontsize=8,
    title_fontsize=9,
    frameon=False
)

plt.tight_layout()

plt.savefig(
    OUTPUT_DIR/"Fig2A_PeptideMotifs.png",
    dpi=600,
    bbox_inches="tight"
)


plt.close()

print("Done.")