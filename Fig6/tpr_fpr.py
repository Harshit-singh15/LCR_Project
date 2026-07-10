import os
import pandas as pd
import sys

# =====================================================
# INPUTS
# =====================================================

metrics_file = sys.argv[1]
confusion_file = sys.argv[2]
output_dir = sys.argv[3]

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

    grouped = (
        df.groupby(
            ["Tool", column_name],
            observed=False,
            as_index=False
        )[["TP", "FP", "FN", "TN"]]
        .sum()
    )

    grouped["TPR"] = (
        grouped["TP"] /
        (grouped["TP"] + grouped["FN"])
    ).fillna(0)

    grouped["FPR"] = (
        grouped["FP"] /
        (grouped["FP"] + grouped["TN"])
    ).fillna(0)

    tpr_df = grouped[
        ["Tool", column_name, "TPR"]
    ].rename(
        columns={
            column_name: "Category"
        }
    )

    fpr_df = grouped[
        ["Tool", column_name, "FPR"]
    ].rename(
        columns={
            column_name: "Category"
        }
    )

    tpr_df.to_csv(
        os.path.join(
            output_dir,
            f"{parameter_name}_tpr.tsv"
        ),
        sep="\t",
        index=False
    )

    fpr_df.to_csv(
        os.path.join(
            output_dir,
            f"{parameter_name}_fpr.tsv"
        ),
        sep="\t",
        index=False
    )

    print(f"Saved {parameter_name}")

print("\nDone.")