import pandas as pd
from pathlib import Path

INPUT_DIR = Path(
    "Fig2/04_metrics/peptide_counts"
)

motif_counts = {}

for file in INPUT_DIR.glob(
    "*_peptide_counts.tsv"
):

    df = pd.read_csv(
        file,
        sep="\t"
    )

    for _, row in df.iterrows():

        motif = row["Best-Peptide"]

        motif_counts[motif] = (
            motif_counts.get(motif,0)
            + row["Count"]
        )

top = (
    pd.Series(motif_counts)
      .sort_values(
          ascending=False
      )
)

print(
    top.head(40)
)