import pandas as pd

df = pd.read_csv(
    "bed_files\\flps_strict_masked_mouse.bed",
    sep="\t"
)

print(df.head(30))