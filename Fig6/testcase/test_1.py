import pandas as pd

df = pd.read_csv(
    r"Fig6\reference_metrics.tsv",
    sep="\t"
)

print(df.columns)

print(df["Length_Bin"].value_counts(dropna=False))

print(df["Count_Bin"].value_counts(dropna=False))

print(df["Entropy_Bin"].value_counts(dropna=False))

print(df["Coverage_Bin"].value_counts(dropna=False))