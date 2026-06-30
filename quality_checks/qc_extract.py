from pathlib import Path
import pandas as pd

INPUT_DIR = Path("outputs_prerequisite/extracted_sequences")

print("\n=== EXTRACTION QC ===\n")

for file in INPUT_DIR.glob("*_lcrs.tsv"):

    df = pd.read_csv(file, sep="\t")

    bad = (
        df["Length"]
        !=
        (df["End"] - df["Start"] + 1)
    ).sum()

    print(
        f"{file.stem:40s} "
        f"Rows={len(df):8d} "
        f"Bad={bad}"
    )