import pandas as pd

df = pd.read_csv(
    "outputs_prerequisite/extracted_sequences/flps_strict_masked_mouse_lcrs.tsv",
    sep="\t"
)

print(df["Length"].describe())

print(df["Length"].value_counts().head(20))