import pandas as pd

# 1. Load the TSV file
# Replace the path with your actual file path if it differs
file_path = "Fig5\\02_metrics\\alcor_mode2_metrics.tsv"
df = pd.read_csv(file_path, sep="\t")

# 2. Display the descriptive statistics for the "Purity" column
print(df.describe())