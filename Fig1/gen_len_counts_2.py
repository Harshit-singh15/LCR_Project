from pathlib import Path
import pandas as pd
import sys
# ==========================================
# INPUT / OUTPUT
# ==========================================

INPUT_DIR = sys.argv[1]
PROTEIN_LENGTH_FILE = sys.argv[2]
LENGTH_DIR = sys.argv[3]
COUNT_DIR = sys.argv[4]

 # LCR counts per protein

Path(LENGTH_DIR).mkdir(parents=True, exist_ok=True)
Path(COUNT_DIR).mkdir(parents=True, exist_ok=True)

# ==========================================
# TOTAL PROTEINS 
# ==========================================
protein_lengths = pd.read_csv(
    PROTEIN_LENGTH_FILE,
    sep="\t"
)

TOTAL_PROTEINS = len(protein_lengths)

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

    df = pd.read_csv(file, sep="\t")

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

        length_results.append([
            label,
            int(count)
        ])

    length_df = pd.DataFrame(
        length_results,
        columns=["Category", "Count"]
    )

    length_output = (
        Path(LENGTH_DIR)
        /
        f"{file.stem.replace('_lcrs','')}_categorized.tsv"
    )

    length_df.to_csv(
        length_output,
        sep="\t",
        index=False
    )

    # =====================================================
    # FIGURE 1C : LCR COUNTS PER PROTEIN
    # =====================================================

    protein_counts = (
        df.groupby("Protein_ID")
          .size()
          .reset_index(name="LCR_Count")
    )

    proteins_with_lcr = len(protein_counts)

    zero_lcr = TOTAL_PROTEINS - proteins_with_lcr

    count_results = []

    # 0 LCRs

    count_results.append([
        "0",
        int(zero_lcr)
    ])

    # 1–5

    count_results.append([
        "1-5",
        int(
            (
                (protein_counts["LCR_Count"] >= 1) &
                (protein_counts["LCR_Count"] <= 5)
            ).sum()
        )
    ])

    # 6–10

    count_results.append([
        "6-10",
        int(
            (
                (protein_counts["LCR_Count"] >= 6) &
                (protein_counts["LCR_Count"] <= 10)
            ).sum()
        )
    ])

    # 11–15

    count_results.append([
        "11-15",
        int(
            (
                (protein_counts["LCR_Count"] >= 11) &
                (protein_counts["LCR_Count"] <= 15)
            ).sum()
        )
    ])

    # 16+

    count_results.append([
        "16+",
        int(
            (
                protein_counts["LCR_Count"] >= 16
            ).sum()
        )
    ])

    count_df = pd.DataFrame(
        count_results,
        columns=["Category", "Count"]
    )

    count_output = (
        Path(COUNT_DIR)
        /
        f"{file.stem.replace('_lcrs','')}_categorized.tsv"
    )

    count_df.to_csv(
        count_output,
        sep="\t",
        index=False
    )

    # =====================================================
    # QC CHECK
    # =====================================================

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

print("\nDone.")