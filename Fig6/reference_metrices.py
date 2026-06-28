import pandas as pd
import numpy as np
import math
from collections import Counter
from Bio import SeqIO

# ===== EDIT THESE =====

reference_bed = r"arabidopsis\dataforFig6\arabidopsis_windows_real.bed"

proteome_fasta = r"arabidopsis\arabidopsis_cleaned.fasta"

output_file = r"arabidopsis\dataforFig6\reference_metrics.tsv"

# ======================


def shannon_entropy(sequence):

    if len(sequence) == 0:
        return 0.0

    counts = Counter(sequence)

    n = len(sequence)

    return -sum(
        (count / n) * math.log2(count / n)
        for count in counts.values()
    )


print("Loading proteome...")

proteins = {}

for record in SeqIO.parse(proteome_fasta, "fasta"):

    proteins[record.id] = str(record.seq)

print(f"Loaded {len(proteins):,} proteins")

print("Loading reference BED...")

df = pd.read_csv(
    reference_bed,
    sep="\t"
)

# Keep only LCR regions

lcr_df = df[
    df["Classification"] == "LCR"
].copy()

print(f"LCR regions: {len(lcr_df):,}")

results = []

protein_ids = set(proteins.keys())

for protein_id in protein_ids:

    sequence = proteins[protein_id]

    protein_length = len(sequence)

    protein_entropy = shannon_entropy(sequence)

    regions = lcr_df[
        lcr_df["Protein_ID"] == protein_id
    ]

    lcr_count = len(regions)

    if lcr_count == 0:

        coverage_percent = 0.0

        entropy_ratio = 0.0

    else:

        lcr_length = (
            regions["End_Position"]
            - regions["Start_Position"]
            + 1
        ).sum()

        coverage_percent = (
            lcr_length
            / protein_length
        ) * 100

        lcr_sequence = []

        for _, row in regions.iterrows():

            start = int(row["Start_Position"])
            end = int(row["End_Position"])

            lcr_sequence.append(
                sequence[start - 1:end]
            )

        lcr_sequence = "".join(
            lcr_sequence
        )

        lcr_entropy = shannon_entropy(
            lcr_sequence
        )

        if protein_entropy > 0:

            entropy_ratio = (
                lcr_entropy
                / protein_entropy
            )

        else:

            entropy_ratio = 0.0

    results.append(
        [
            protein_id,
            protein_length,
            lcr_count,
            coverage_percent,
            entropy_ratio
        ]
    )

metrics = pd.DataFrame(
    results,
    columns=[
        "Protein_ID",
        "Protein_Length",
        "LCR_Count",
        "Coverage_Percent",
        "Entropy_Ratio"
    ]
)

# --------------------
# Length bins (1-10)
# --------------------

metrics["Length_Bin"] = pd.qcut(
    metrics["Protein_Length"],
    q=10,
    labels=[str(i) for i in range(1,11)],
    duplicates="drop"
)

# --------------------
# LCR count bins
# --------------------

def count_bin(x):

    if x >= 6:
        return "6"

    return str(int(x))

metrics["Count_Bin"] = metrics[
    "LCR_Count"
].apply(count_bin)

# --------------------
# Coverage bins
# --------------------

def coverage_bin(x):

    if x <= 5:
        return "5"

    elif x <= 10:
        return "10"

    elif x <= 15:
        return "15"

    elif x <= 20:
        return "20"

    else:
        return ">20"

metrics["Coverage_Bin"] = metrics[
    "Coverage_Percent"
].apply(coverage_bin)

# --------------------
# Entropy bins
# --------------------

def entropy_bin(x):

    if x <= 0.2:
        return "0.2"

    elif x <= 0.4:
        return "0.4"

    elif x <= 0.6:
        return "0.6"

    elif x <= 0.8:
        return "0.8"

    else:
        return "1.0"

metrics["Entropy_Bin"] = metrics[
    "Entropy_Ratio"
].apply(entropy_bin)

metrics.to_csv(
    output_file,
    sep="\t",
    index=False
)

print(f"\nSaved: {output_file}")

print(metrics.head())