import pandas as pd
from pathlib import Path

INPUT_DIR = Path(
    "Fig2/04_metrics/filtered"
)

OUTPUT_DIR = Path(
    "Fig2/04_metrics/peptide_counts"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

for file in INPUT_DIR.glob(
    "*_filtered.tsv"
):

    df = pd.read_csv(
        file,
        sep="\t"
    )

    counts = (
        df["Best-Peptide"]
        .value_counts()
        .reset_index()
    )

    counts.columns = [
        "Best-Peptide",
        "Count"
    ]

    counts["Proportion"] = (
        counts["Count"] /
        counts["Count"].sum()
    )

    out_file = (
        OUTPUT_DIR /
        file.name.replace(
            "_filtered.tsv",
            "_peptide_counts.tsv"
        )
    )

    counts.to_csv(
        out_file,
        sep="\t",
        index=False
    )

    print(
        f"Saved: {out_file}"
    )