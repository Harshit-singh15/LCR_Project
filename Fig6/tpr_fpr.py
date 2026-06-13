import os
import pandas as pd

# =====================================================
# INPUTS
# =====================================================

metrics_file = r"Fig6\reference_metrics.tsv"

confusion_file = r"Fig6\protein_confusion.tsv"

output_dir = r"Fig6\plot_tables"

# =====================================================

os.makedirs(output_dir, exist_ok=True)

print("Loading files...")

metrics = pd.read_csv(
    metrics_file,
    sep="\t"
)

confusion = pd.read_csv(
    confusion_file,
    sep="\t"
)

print("Merging...")

df = confusion.merge(
    metrics,
    on="Protein_ID",
    how="left"
)

print(df.shape)

# =====================================================
# PARAMETERS
# =====================================================

parameter_map = {
    "gene_length": "Length_Bin",
    "lcr_count": "Count_Bin",
    "coverage": "Coverage_Bin",
    "entropy_ratio": "Entropy_Bin"
}

# =====================================================
# BUILD TABLES
# =====================================================

for parameter_name, column_name in parameter_map.items():

    print(f"\nProcessing {parameter_name}")

    rows_tpr = []
    rows_fpr = []

    grouped = df.groupby(
        ["Tool", column_name],
        observed=False
    )

    for (tool, category), group in grouped:

        tp = group["TP"].sum()
        fp = group["FP"].sum()
        fn = group["FN"].sum()
        tn = group["TN"].sum()

        if (tp + fn) > 0:
            tpr = tp / (tp + fn)
        else:
            tpr = 0

        if (fp + tn) > 0:
            fpr = fp / (fp + tn)
        else:
            fpr = 0

        rows_tpr.append(
            [
                tool,
                category,
                tpr
            ]
        )

        rows_fpr.append(
            [
                tool,
                category,
                fpr
            ]
        )

    tpr_df = pd.DataFrame(
        rows_tpr,
        columns=[
            "Tool",
            "Category",
            "TPR"
        ]
    )

    fpr_df = pd.DataFrame(
        rows_fpr,
        columns=[
            "Tool",
            "Category",
            "FPR"
        ]
    )

    tpr_file = os.path.join(
        output_dir,
        f"{parameter_name}_tpr.tsv"
    )

    fpr_file = os.path.join(
        output_dir,
        f"{parameter_name}_fpr.tsv"
    )

    tpr_df.to_csv(
        tpr_file,
        sep="\t",
        index=False
    )

    fpr_df.to_csv(
        fpr_file,
        sep="\t",
        index=False
    )

    print(
        f"Saved {parameter_name}"
    )

print("\nDone.")