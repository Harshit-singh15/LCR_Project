import pandas as pd
import glob

rows = []

for f in glob.glob("Fig2/02_consensus_beds/*.bed"):

    n = int(
        f.split("_")[-1].replace(".bed","")
    )

    count = sum(1 for _ in open(f))

    rows.append([n,count])

summary = pd.DataFrame(
    rows,
    columns=["Consensus","Count"]
)

summary.sort_values(
    "Consensus"
).to_csv(
    "Fig2/consensus_counts.tsv",
    sep="\t",
    index=False
)