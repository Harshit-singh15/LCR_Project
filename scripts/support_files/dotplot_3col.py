import os
import pandas as pd

# 1. Configuration - Change these filenames to match yours
INPUT_FILE = "bed_files\\Dotplot.bed"  # Your existing 4-column BED file
OUTPUT_FILE = "new_three_column.bed"  # The new file that will be created


def extract_bed_columns(infile, outfile):
    if not os.path.exists(infile):
        print(f"Error: The input file '{infile}' does not exist.")
        return

    print(f"Reading {infile}...")

    # Read the file (headerless, tab or space-separated)
    # usecols=[0, 1, 2] ensures we only load the first 3 columns into memory
    df = pd.read_csv(
        infile, sep=r"\s+", header=None, usecols=[0, 1, 2], engine="python"
    )

    # Save to a completely new file
    # sep='\t', header=False, index=False ensures valid BED format layout
    df.to_csv(outfile, sep="\t", header=False, index=False)

    print(f"Done! Created a new 3-column BED file at: '{outfile}'")


if __name__ == "__main__":
    extract_bed_columns(INPUT_FILE, OUTPUT_FILE)