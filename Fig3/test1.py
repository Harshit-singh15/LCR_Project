import pandas as pd

# 1. Load the TSV file
# Replace the path with your actual file path if it differs
file_path = "Fig3\\02_metrics\\Dotplot_purity.tsv"
df = pd.read_csv(file_path, sep="\t")

# 2. Display the descriptive statistics for the "Purity" column
print(df["Purity"].describe())