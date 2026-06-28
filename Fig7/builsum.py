import pandas as pd
from pathlib import Path

input_dir = Path(r"arabidopsis\dataforFig6\plot_tables")
output_dir = Path(r"arabidopsis\dataforFig7")

output_dir.mkdir(exist_ok=True)

datasets = {
    "gene_length":
        ("gene_length_tpr.tsv",
         "gene_length_fpr.tsv"),

    "lcr_count":
        ("lcr_count_tpr.tsv",
         "lcr_count_fpr.tsv"),

    "coverage":
        ("coverage_tpr.tsv",
         "coverage_fpr.tsv"),

    "entropy_ratio":
        ("entropy_ratio_tpr.tsv",
         "entropy_ratio_fpr.tsv")
}

for name, files in datasets.items():

    tpr = pd.read_csv(
        input_dir / files[0],
        sep="\t"
    )

    fpr = pd.read_csv(
        input_dir / files[1],
        sep="\t"
    )

    mean_tpr = (
        tpr.groupby("Tool")["TPR"]
        .mean()
        .reset_index()
    )

    mean_fpr = (
        fpr.groupby("Tool")["FPR"]
        .mean()
        .reset_index()
    )

    merged = mean_tpr.merge(
        mean_fpr,
        on="Tool"
    )

    merged.to_csv(
        output_dir / f"{name}_summary.tsv",
        sep="\t",
        index=False
    )

print("Done")