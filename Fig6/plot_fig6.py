import os
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

# =====================================================
# INPUT AND OUTPUT DIRECTORIES
# =====================================================
input_dir = Path(r"zebrafish\dataforFig6\plot_tables")
output_dir = Path(r"zebrafish\Fig_outputs\Fig6")

# Ensure the output directory exists so saving doesn't fail
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
def plot_metric(ax, df, value_col, title, xlab, order):
    # Drop rows missing critical categorization info to prevent crash
    df = df.dropna(subset=["Category", value_col]).copy()

    tools = sorted(df["Tool"].unique())

    for tool in tools:
        temp = df[df["Tool"] == tool].copy()
        
        # Safely normalize float-like integers (e.g., 5.0 -> 5) before converting to str
        temp["Category"] = temp["Category"].apply(
            lambda x: str(int(x)) if isinstance(x, (int, float)) and x == int(x) else str(x)
        )

        # Enforce categorical sequence strictly mapping to strings
        categories_str = [str(x) for x in order]
        temp["Category"] = pd.Categorical(
            temp["Category"],
            categories=categories_str,
            ordered=True
        )

        # Drop values that fall completely out of your predefined order bucket list
        temp = temp.dropna(subset=["Category"])
        temp = temp.sort_values("Category")

        # Explicitly conversion to string for Matplotlib line categorical mapping
        x_data = temp["Category"].dt.categories if hasattr(temp["Category"], "dt") else temp["Category"].astype(str)

        ax.plot(
            x_data,
            temp[value_col],
            marker="o",
            linewidth=1.5,
            markersize=4,
            label=tool
        )

    ax.set_title(title, fontsize=11)
    ax.set_xlabel(xlab)
    ax.set_ylabel(value_col)
    ax.set_ylim(0, 1)
    ax.set_yticks(
        np.arange(
        0,
        1.00,
        0.25
        )
    )
    ax.tick_params(axis="x", rotation=45)
    ax.grid(alpha=0.3)

# =====================================================
# GENERATE PANELS
# =====================================================

# --- PANEL A ---
fig, axes = plt.subplots(1, 2, figsize=(14, 8), constrained_layout=True)
plot_metric(axes[0], gene_tpr, "TPR", "Gene Length Categories", "Length Bin", gene_order)
plot_metric(axes[1], gene_fpr, "FPR", "Gene Length Categories", "Length Bin", gene_order)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=8)
plt.savefig(output_dir / "Fig6A_GeneLength.png", dpi=300, bbox_inches="tight")
plt.close()

# --- PANEL B ---
fig, axes = plt.subplots(1, 2, figsize=(14, 8), constrained_layout=True)
plot_metric(axes[0], count_tpr, "TPR", "LCR Number Per Gene", "LCR Count", count_order)
plot_metric(axes[1], count_fpr, "FPR", "LCR Number Per Gene", "LCR Count", count_order)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=8)
plt.savefig(output_dir / "Fig6B_LCRCount.png", dpi=300, bbox_inches="tight")
plt.close()

# --- PANEL C ---
fig, axes = plt.subplots(1, 2, figsize=(14, 8), constrained_layout=True)
plot_metric(axes[0], coverage_tpr, "TPR", "LCR Coverage (%)", "Coverage", coverage_order)
plot_metric(axes[1], coverage_fpr, "FPR", "LCR Coverage (%)", "Coverage", coverage_order)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=8)
plt.savefig(output_dir / "Fig6C_Coverage.png", dpi=300, bbox_inches="tight")
plt.close()

# --- PANEL D ---
fig, axes = plt.subplots(1, 2, figsize=(14, 8), constrained_layout=True)
plot_metric(axes[0], entropy_tpr, "TPR", "LCR : Gene Entropy Ratio", "Entropy Ratio", entropy_order)
plot_metric(axes[1], entropy_fpr, "FPR", "LCR : Gene Entropy Ratio", "Entropy Ratio", entropy_order)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=8)
plt.savefig(output_dir / "Fig6D_EntropyRatio.png", dpi=300, bbox_inches="tight")
plt.close()

print("All panels successfully constructed and saved under Fig6/plots/")