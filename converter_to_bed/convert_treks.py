from pathlib import Path
import sys
import pandas as pd

print("[INFO] Running convert_trekes")

# 1. PASTE YOUR EXACT PATHS HERE
# Use absolute paths (e.g., r"C:\path\to\file.tsv") or relative paths.
input_file = Path(sys.argv[1])  # Get the input file path from command line argument
output_file = Path(sys.argv[2])  # Get the output file path from command line argument

# Ensure the parent directory for the output file exists
output_file.parent.mkdir(parents=True, exist_ok=True)

# Validate that the input file exists before running
if not input_file.is_file():
    raise FileNotFoundError(f"Input file not found: {input_file}")

print(f"[INFO] Processing {input_file.name}")

# Read the single TSV file
df = pd.read_csv(input_file, sep="\t")

# Extract UniProt accession
def get_accession(seqid):
    parts = str(seqid).split("|")
    if len(parts) >= 3:
        return parts[1]
    return seqid

df["Protein_ID"] = df["seqid"].apply(get_accession)

# Write out to the specific BED file path
df[["Protein_ID", "start", "end"]].rename(
    columns={
        "start": "Start",
        "end": "End"
    }
).to_csv(
    output_file,
    sep="\t",
    index=False,
    header=False
)

print(f"[SUCCESS] {input_file.name} converted")
print("T-REKS BED file written:", output_file)
print("Regions:", len(df))
print("[SUCCESS] Converter completed")