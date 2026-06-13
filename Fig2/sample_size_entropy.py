import pandas as pd

e = pd.read_csv(
    "Fig2/04_metrics/entropy.tsv",
    sep="\t"
)

print(
    e.groupby("Consensus")
     .size()
)