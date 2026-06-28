import os
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# =====================================================
# INPUT AND OUTPUT DIRECTORIES
# =====================================================
input_dir = Path(r"fruitfly\dataforFig6\plot_tables")
output_dir = Path(r"fruitfly\Fig_outputs\Fig6")

# Ensure the output directory exists
output_dir.mkdir(parents=True, exist_ok=True)

# =====================================================
# LOAD FILES
# =====================================================
gene_tpr = pd.read_csv(input_dir / "gene_length_tpr.tsv", sep="\t")
gene_fpr = pd.read_csv(input_dir / "gene_length_fpr.tsv", sep="\t")

count_tpr = pd.read_csv(input_dir / "lcr_count_tpr.tsv", sep="\t")
count_fpr = pd.read_csv(input_dir / "lcr_count_fpr.tsv", sep="\t")

coverage_tpr = pd.read_csv(input_dir / "coverage_tpr.tsv", sep="\t")
coverage_fpr = pd.read_csv(input_dir / "coverage_fpr.tsv", sep="\t")

entropy_tpr = pd.read_csv(input_dir / "entropy_ratio_tpr.tsv", sep="\t")
entropy_fpr = pd.read_csv(input_dir / "entropy_ratio_fpr.tsv", sep="\t")

# =====================================================
# CATEGORY ORDERING
# =====================================================
gene_order = [str(i) for i in range(1, 11)]
count_order = ["1", "2", "3", "4", "5", "6"]
coverage_order = ["5", "10", "15", "20", ">20"]
entropy_order = ["0.2", "0.4", "0.6", "0.8", "1.0"]

# =====================================================
# PLOTTING FUNCTION
# =====================================================
def plot_metric(ax, df, value_col, title, xlab, order, panel_label):
    # Data Cleanup to avoid data type float errors
    df = df.dropna(subset=["Category", value_col]).copy()
    tools = sorted(df["Tool"].unique())

    for tool in tools:
        temp = df[df["Tool"] == tool].copy()
        
        # Safely convert category rows into cleanly normalized string indexes
        temp["Category"] = temp["Category"].apply(
            lambda x: str(int(x)) if isinstance(x, (int, float)) and x == int(x) else str(x)
        )

        categories_str = [str(x) for x in order]
        temp["Category"] = pd.Categorical(
            temp["Category"],
            categories=categories_str,
            ordered=True
        )

        temp = temp.dropna(subset=["Category"])
        temp = temp.sort_values("Category")

        # Fallback categorical check to satisfy Matplotlib requirements
        x_data = temp["Category"].astype(str)

        ax.plot(
            x_data,
            temp[value_col],
            marker="o",
            linewidth=1.5,
            markersize=4,
            label=tool
        )

    # Stylistic setups 
    ax.set_title(title, fontsize=12, fontweight="bold", loc="center")
    ax.set_xlabel(xlab, fontsize=10)
    ax.set_ylabel(value_col, fontsize=10)
    ax.set_ylim(0, 1.05)
    ax.tick_params(axis="x", rotation=45, labelsize=9)
    ax.tick_params(axis="y", labelsize=9)
    ax.grid(alpha=0.3)
    
    # Add Panel Label Indicator tag to the top-left corner of the sub-block
    ax.text(-0.15, 1.1, panel_label, transform=ax.transAxes, 
            fontsize=16, fontweight="bold", va="top", ha="right")

# =====================================================
# MASTER PLOT SETUP (4 Rows x 2 Columns)
# =====================================================
# gridspec_kw adds strategic bottom padding workspace space allocated for the global shared legend box.
fig, axes = plt.subplots(
    nrows=4, 
    ncols=2, 
    figsize=(14, 32), 
    gridspec_kw={'bottom': 0.06},
    layout="constrained"
)

# Row 0: Panel A (Gene Length)
plot_metric(axes[0, 0], gene_tpr, "TPR", "Gene Length Categories", "Length Bin", gene_order, "A")
plot_metric(axes[0, 1], gene_fpr, "FPR", "Gene Length Categories", "Length Bin", gene_order, "")

# Row 1: Panel B (LCR Count)
plot_metric(axes[1, 0], count_tpr, "TPR", "LCR Number Per Gene", "LCR Count", count_order, "B")
plot_metric(axes[1, 1], count_fpr, "FPR", "LCR Number Per Gene", "LCR Count", count_order, "")

# Row 2: Panel C (LCR Coverage)
plot_metric(axes[2, 0], coverage_tpr, "TPR", "LCR Coverage (%)", "Coverage", coverage_order, "C")
plot_metric(axes[2, 1], coverage_fpr, "FPR", "LCR Coverage (%)", "Coverage", coverage_order, "")

# Row 3: Panel D (Entropy Ratio)
plot_metric(axes[3, 0], entropy_tpr, "TPR", "LCR : Gene Entropy Ratio", "Entropy Ratio", entropy_order, "D")
plot_metric(axes[3, 1], entropy_fpr, "FPR", "LCR : Gene Entropy Ratio", "Entropy Ratio", entropy_order, "")

# =====================================================
# COMPILING SINGLE SHARED GLOBAL LEGEND
# =====================================================
# Extract legend labels from one of the plots
handles, labels = axes[0, 0].get_legend_handles_labels()

# Draw one legend down across the center at the bottom of the canvas
fig.legend(
    handles,
    labels,
    loc="lower center",
    ncol=min(len(labels), 6),
    fontsize=11,
    frameon=True,
    facecolor="white",
    edgecolor="lightgray"
)

# =====================================================
# SAVE FILES
# =====================================================
plt.savefig(output_dir / "Figure6_Combined.png", dpi=300, bbox_inches="tight")

plt.close()
print("Success! Master 4x2 grid saved as 'Figure6_Combined.png")