from pathlib import Path
from collections import Counter
import pandas as pd

# =====================================================
# CONFIG
# =====================================================

INPUT_DIR = "celegans\\outputs_prerequisite\\extracted_sequences"
OUTPUT_DIR = "celegans\\outputs_prerequisite\\Amino_acid"

Path(OUTPUT_DIR).mkdir(
    parents=True,
    exist_ok=True
)

# Standard amino acids

AA_ORDER = list(
    "ACDEFGHIKLMNPQRSTVWY"
)

# =====================================================
# PROCESS FILES
# =====================================================

for file in Path(INPUT_DIR).glob("*_lcrs.tsv"):

    print(f"\nProcessing {file.name}")

    df = pd.read_csv(
        file,
        sep="\t"
    )

    aa_counts = Counter()

    total_residues = 0

    for seq in df["Sequence"].fillna(""):

        seq = str(seq)

        aa_counts.update(seq)

        total_residues += len(seq)

    records = []

    for aa in AA_ORDER:

        count = aa_counts.get(
            aa,
            0
        )

        proportion = (
            count /
            total_residues
            if total_residues > 0
            else 0
        )

        records.append(
            [
                aa,
                count,
                proportion
            ]
        )

    out_df = pd.DataFrame(
        records,
        columns=[
            "Character",
            "Count",
            "Proportion"
        ]
    )

    print(
        f"Total residues: {total_residues}"
    )

    print(
        f"Sum proportions: "
        f"{out_df['Proportion'].sum():.6f}"
    )

    out_file = (
        Path(OUTPUT_DIR)
        /
        f"{file.stem.replace('_lcrs','')}_aa_count"
    )

    out_df.to_csv(
        out_file,
        sep="\t",
        index=False
    )

print("\nAA composition complete.")