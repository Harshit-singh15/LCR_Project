import pandas as pd
from pathlib import Path

folder = Path(
    "Fig2/04_metrics/peptide_counts"
)

all_counts = {}

for file in folder.glob("*_peptide_counts.tsv"):

    df = pd.read_csv(
        file,
        sep="\t"
    )

    for _, row in df.iterrows():

        motif = row["Best-Peptide"]

        all_counts[motif] = (
            all_counts.get(motif, 0)
            + row["Count"]
        )

top_motifs = (
    pd.Series(all_counts)
      .sort_values(
          ascending=False
      )
      .head(15)
)

print(top_motifs)