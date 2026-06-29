from pathlib import Path
from collections import Counter
import pandas as pd
import math

INPUT_DIR = r"ecoli\dataforFig1\extracted_sequences"
OUTPUT_DIR = r"ecoli\dataforFig1\ShanonEntropy"

Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

def shannon_entropy(seq):

    if not seq:
        return 0

    counts = Counter(seq)

    entropy = 0.0

    for count in counts.values():

        p = count / len(seq)

        entropy -= p * math.log2(p)

    return entropy

for file in Path(INPUT_DIR).glob("*_lcrs.tsv"):

    print(f"\nProcessing {file.name}")

    df = pd.read_csv(file, sep="\t")

    entropies = [
        shannon_entropy(str(seq))
        for seq in df["Sequence"]
    ]

    out_df = pd.DataFrame({
        "Shannon_Entropy": entropies
    })

    print(
        f"Mean = {out_df['Shannon_Entropy'].mean():.3f}"
    )

    print(
        f"Max = {out_df['Shannon_Entropy'].max():.3f}"
    )

    out_file = (
        Path(OUTPUT_DIR)
        /
        f"{file.stem.replace('_lcrs','')}_SNS"
    )

    out_df.to_csv(
        out_file,
        sep="\t",
        index=False
    )

print("\nDone.")