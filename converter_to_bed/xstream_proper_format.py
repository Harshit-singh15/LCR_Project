import pandas as pd

input_file = "bed_files\\xstream.bed"
output_file = "xstream_fixed.bed"

df = pd.read_csv(input_file, sep="\t")

# Keep only the UniProt identifier part
df.iloc[:, 0] = (
    df.iloc[:, 0]
    .astype(str)
    .str.split()
    .str[0]
)

# Keep only first 3 columns
df = df.iloc[:, :3]

df.columns = ["Protein_ID", "Start", "End"]

df.to_csv(
    output_file,
    sep="\t",
    index=False
)

print(f"Saved: {output_file}")