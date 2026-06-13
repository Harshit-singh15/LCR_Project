from pathlib import Path
import pandas as pd

INPUT_DIR = Path("outputs_prerequisite/Diversity")

print("\n=== ENTROPY QC ===\n")

for file in INPUT_DIR.glob("*_SNS"):

    df = pd.read_csv(file, sep="\t")

    mean_e = df["Shannon_Entropy"].mean()
    max_e = df["Shannon_Entropy"].max()

    print(
        f"{file.stem:40s} "
        f"Mean={mean_e:.3f} "
        f"Max={max_e:.3f}"
    )