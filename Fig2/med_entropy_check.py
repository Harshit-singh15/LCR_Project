import pandas as pd
import matplotlib.pyplot as plt

e = pd.read_csv(
    "Fig2/04_metrics/entropy.tsv",
    sep="\t"
)

med = (
    e.groupby("Consensus")["Entropy"]
    .median()
)

plt.figure(figsize=(8,5))

plt.plot(
    med.index,
    med.values,
    marker="o"
)

plt.xlabel("Consensus")
plt.ylabel("Median Entropy")

plt.show()