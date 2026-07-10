from pathlib import Path
import pandas as pd
import gc
import sys

# ==========================================
# INPUT / OUTPUT
# ==========================================

INPUT_DIR = sys.argv[1]
PROTEIN_LENGTH_FILE = sys.argv[2]
LENGTH_DIR = sys.argv[3]
COUNT_DIR = sys.argv[4]

Path(LENGTH_DIR).mkdir(parents=True, exist_ok=True)
Path(COUNT_DIR).mkdir(parents=True, exist_ok=True)

# ==========================================
# TOTAL PROTEINS
# ==========================================

TOTAL_PROTEINS = sum(1 for _ in open(PROTEIN_LENGTH_FILE)) - 1

# ==========================================
# LENGTH BINS (Figure 1A)
# ==========================================

length_bins = [
    ("0-10",    0,   10),
    ("10-20",   10,  20),
    ("20-50",   20,  50),
    ("50-100",  50, 100),
    ("100-200", 100, 200),
    ("200+",    200, float("inf"))
]

# ==========================================
# PROCESS EACH TOOL
# ==========================================

for file in Path(INPUT_DIR).glob("*_lcrs.tsv"):

    print(f"\nProcessing {file.name}")

    df = pd.read_csv(
        file,
        sep="\t",
        usecols=["Protein_ID", "Length"],
        dtype={
            "Protein_ID": "string",
            "Length": "int32"
        }
    )

    # =====================================================
    # FIGURE 1A : LCR LENGTH DISTRIBUTION
    # =====================================================

    length_results = []

    for label, low, high in length_bins:

        if high == float("inf"):
            count = (df["Length"] > low).sum()
        else:
            count = (
                (df["Length"] > low) &
                (df["Length"] <= high)
            ).sum()

        length_results.append((label, int(count)))

    length_output = (
        Path(LENGTH_DIR)
        /
        f"{file.stem.replace('_lcrs','')}_categorized.tsv"
    )

    pd.DataFrame(
        length_results,
        columns=["Category", "Count"]
    ).to_csv(
        length_output,
        sep="\t",
        index=False
    )

    # =====================================================
    # FIGURE 1C : LCR COUNTS PER PROTEIN
    # =====================================================

    protein_counts = df.groupby(
        "Protein_ID",
        sort=False
    ).size()

    proteins_with_lcr = len(protein_counts)

    zero_lcr = TOTAL_PROTEINS - proteins_with_lcr

    count_results = [

        (
            "0",
            int(zero_lcr)
        ),

        (
            "1-5",
            int(
                (
                    (protein_counts >= 1) &
                    (protein_counts <= 5)
                ).sum()
            )
        ),

        (
            "6-10",
            int(
                (
                    (protein_counts >= 6) &
                    (protein_counts <= 10)
                ).sum()
            )
        ),

        (
            "11-15",
            int(
                (
                    (protein_counts >= 11) &
                    (protein_counts <= 15)
                ).sum()
            )
        ),

        (
            "16+",
            int(
                (
                    protein_counts >= 16
                ).sum()
            )
        )

    ]

    count_output = (
        Path(COUNT_DIR)
        /
        f"{file.stem.replace('_lcrs','')}_categorized.tsv"
    )

    count_df = pd.DataFrame(
        count_results,
        columns=["Category", "Count"]
    )

    count_df.to_csv(
        count_output,
        sep="\t",
        index=False
    )

    total_check = count_df["Count"].sum()

    print(
        f"Length categories saved -> {length_output.name}"
    )

    print(
        f"Count categories saved  -> {count_output.name}"
    )

    print(
        f"Protein total check     -> {total_check}"
    )

    del df
    del protein_counts
    del count_df
    gc.collect()

print("\nDone.")