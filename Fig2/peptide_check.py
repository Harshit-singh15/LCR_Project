import pandas as pd
from pathlib import Path

folder = Path(
    "Fig2/04_metrics/peptide_counts"
)

for file in sorted(
    folder.glob("*_peptide_counts.tsv")
):

    df = pd.read_csv(
        file,
        sep="\t"
    )

    top_prop = df["Proportion"].max()

    print(
        file.stem,
        round(top_prop, 3)
    )