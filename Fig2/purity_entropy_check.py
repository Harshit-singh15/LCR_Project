import pandas as pd

e = pd.read_csv(
    "Fig2/04_metrics/entropy.tsv",
    sep="\t"
)

p = pd.read_csv(
    "Fig2/04_metrics/purity.tsv",
    sep="\t"
)

print(
    e.groupby("Consensus")["Entropy"]
    .median()
)

print()

print(
    p.groupby("Consensus")["Purity"]
    .median()
)