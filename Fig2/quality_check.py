import pandas as pd

e = pd.read_csv(
    "Fig2/04_metrics/entropy.tsv",
    sep="\t"
)

p = pd.read_csv(
    "Fig2/04_metrics/purity.tsv",
    sep="\t"
)

print(e["Entropy"].describe())
print(p["Purity"].describe())