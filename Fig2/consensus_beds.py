import pandas as pd
import os

multi = pd.read_csv(
    "Fig2\\multiinter.tsv",
    sep="\t",
    header=None
)

os.makedirs(
    "Fig2/02_consensus_beds",
    exist_ok=True
)

for n in range(1,14):

    subset = multi[multi[3] == n]

    out = subset[[0,1,2]]

    outfile = (
        f"Fig2/02_consensus_beds/"
        f"consensus_{n}.bed"
    )

    out.to_csv(
        outfile,
        sep="\t",
        header=False,
        index=False
    )

    print(
        f"consensus_{n}:",
        len(out)
    )