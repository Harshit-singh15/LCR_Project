from Bio import SeqIO
from collections import Counter
from pathlib import Path
import pandas as pd

INPUT_DIR = Path("Fig5/01_fastas")
OUTPUT_DIR = Path("Fig5/02_metrics")

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

# ==========================================
# Mutation calculation
# ==========================================

def mutations_to_repeat(seq, k):

    n = len(seq)

    if n < k:
        return n

    total_mutations = 0

    for pos in range(k):

        chars = []

        i = pos

        while i < n:

            chars.append(seq[i])

            i += k

        counts = Counter(chars)

        best = max(counts.values())

        total_mutations += (
            len(chars) - best
        )

    return total_mutations


# ==========================================
# Process each method
# ==========================================

for fasta in INPUT_DIR.glob("*.fa"):

    rows = []

    method = fasta.stem

    print(f"Processing {method}")

    for record in SeqIO.parse(
        fasta,
        "fasta"
    ):

        seq = str(record.seq)

        if len(seq) < 2:
            continue

        # ----------------------
        # Purity
        # ----------------------

        counts = Counter(seq)

        most_freq = max(
            counts.values()
        )

        purity = (
            most_freq /
            len(seq)
        ) * 100

        # ----------------------
        # Mutation %
        # ----------------------

        mono_mut = mutations_to_repeat(
            seq,
            1
        )

        di_mut = mutations_to_repeat(
            seq,
            2
        )

        tri_mut = mutations_to_repeat(
            seq,
            3
        )

        best_mut = min(
            mono_mut,
            di_mut,
            tri_mut
        )

        mutation_pct = (
            best_mut /
            len(seq)
        ) * 100

        rows.append([
            mutation_pct,
            purity
        ])

    df = pd.DataFrame(
        rows,
        columns=[
            "Mutation_Percent",
            "Most_Frequent_AA_Percent"
        ]
    )

    outfile = (
        OUTPUT_DIR /
        f"{method}_metrics.tsv"
    )

    df.to_csv(
        outfile,
        sep="\t",
        index=False
    )

    print(
        f"Saved: {outfile}"
    )

print("Done")