from pathlib import Path
import pandas as pd

INPUT_DIR = Path("outputs_prerequisite/Amino_acid")

print("\n=== AA COMPOSITION QC ===\n")

for file in INPUT_DIR.glob("*_aa_count"):

    df = pd.read_csv(file, sep="\t")

    total = df["Proportion"].sum()

    print(
        f"{file.stem:40s} "
        f"{total:.6f}"
    )