import pandas as pd
import math
from collections import Counter
from Bio import SeqIO
import sys
from pathlib import Path

# =====================================================
# INPUTS
# =====================================================

reference_bed = Path(sys.argv[1])      # Reference BED file
proteome_fasta = Path(sys.argv[2])     # Proteome FASTA
output_file = Path(sys.argv[3])        # Output metrics file

output_file.parent.mkdir(
    parents=True,
    exist_ok=True
)

# =====================================================
# FUNCTIONS
# =====================================================

def shannon_entropy(sequence):

    if not sequence:
        return 0.0

    counts = Counter(sequence)

    n = len(sequence)

    return -sum(
        (count / n) * math.log2(count / n)
        for count in counts.values()
    )

# =====================================================
# LOAD PROTEOME
# =====================================================

print("Loading proteome...")

proteins = {
    record.id: str(record.seq)
    for record in SeqIO.parse(proteome_fasta, "fasta")
}

print(f"Loaded {len(proteins):,} proteins")

# =====================================================
# LOAD REFERENCE BED
# =====================================================

print("Loading reference BED...")

df = pd.read_csv(
    reference_bed,
    sep="\t"
)

lcr_df = df[
    df["Classification"] == "LCR"
].copy()

print(f"LCR regions: {len(lcr_df):,}")

# =====================================================
# GROUP ONCE (Huge CPU improvement)
# =====================================================

grouped_regions = {
    protein: group
    for protein, group in lcr_df.groupby("Protein_ID")
}

# =====================================================
# PROCESS PROTEINS
# =====================================================

results = []

for protein_id, sequence in proteins.items():

    protein_length = len(sequence)

    protein_entropy = shannon_entropy(sequence)

    regions = grouped_regions.get(protein_id)

    if regions is None:

        lcr_count = 0
        coverage_percent = 0.0
        entropy_ratio = 0.0

    else:

        lcr_count = len(regions)

        lcr_length = (
            regions["End_Position"]
            - regions["Start_Position"]
            + 1
        ).sum()

        coverage_percent = (
            lcr_length / protein_length
        ) * 100

        lcr_sequence = []

        for row in regions.itertuples(index=False):

            lcr_sequence.append(
                sequence[
                    row.Start_Position - 1:
                    row.End_Position
                ]
            )

        lcr_sequence = "".join(lcr_sequence)

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

# =====================================================
# DATAFRAME
# =====================================================

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

# =====================================================
# LENGTH BINS
# =====================================================

metrics["Length_Bin"] = pd.qcut(
    metrics["Protein_Length"],
    q=10,
    labels=[str(i) for i in range(1, 11)],
    duplicates="drop"
)

# =====================================================
# LCR COUNT BINS
# =====================================================

def count_bin(x):

    if x >= 6:
        return "6"

    return str(int(x))

metrics["Count_Bin"] = metrics[
    "LCR_Count"
].apply(count_bin)

# =====================================================
# COVERAGE BINS
# =====================================================

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

# =====================================================
# ENTROPY BINS
# =====================================================

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

# =====================================================
# SAVE
# =====================================================

metrics.to_csv(
    output_file,
    sep="\t",
    index=False
)

print(f"\nSaved: {output_file}")

print(metrics.head())