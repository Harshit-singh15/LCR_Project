import pandas as pd

df = pd.read_csv(
    "treks_clustalw.bed",
    sep="\t",
    header=None
)

df = df.dropna(subset=[1, 2])

df[1] = df[1].astype(int)
df[2] = df[2].astype(int)

df.to_csv(
    "clustalw_fixed.bed",
    sep="\t",
    header=False,
    index=False
)