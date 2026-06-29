from pathlib import Path
from collections import Counter
from Bio import SeqIO
import pandas as pd
import math
import re

# ======================================================
# Configuration
# ======================================================

INPUT_DIR = Path(r"ecoli\dataforFig3\01_fastas")
OUTPUT_DIR = Path(r"ecoli\dataforFig5\02_metrics")

# Maximum repeat unit length to test
MAX_K = 6

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

# ======================================================
# Header parser
# ======================================================

header_pattern = re.compile(
    r"(.+):(\d+)-(\d+)"
)

# ======================================================
# Shannon entropy
# ======================================================

def shannon_entropy(seq):

    counts = Counter(seq)

    length = len(seq)

    entropy = 0

    for count in counts.values():

        p = count / length

        entropy -= p * math.log2(p)

    return entropy

# ======================================================
# Mutation calculation
# ======================================================

def mutations_to_repeat(seq, k):

    """
    Minimum mutations required to convert
    the sequence into a perfect k-periodic repeat.
    """

    n = len(seq)

    if n < k:
        return n

    total_mut = 0

    for pos in range(k):

        column = []

        i = pos

        while i < n:

            column.append(seq[i])

            i += k

        counts = Counter(column)

        total_mut += (
            len(column)
            - max(counts.values())
        )

    return total_mut

# ======================================================
# Process each FASTA
# ======================================================

for fasta in sorted(INPUT_DIR.glob("*.fa")):

    method = fasta.stem

    print(f"Processing {method}")

    rows = []

    for record in SeqIO.parse(fasta, "fasta"):

        seq = str(record.seq).upper()

        if len(seq) == 0:
            continue

        # ------------------------------------------
        # Parse header
        # ------------------------------------------

        match = header_pattern.match(record.id)

        if match:

            protein, start, end = match.groups()

        else:

            protein = record.id

            start = ""

            end = ""

        # ------------------------------------------
        # Length
        # ------------------------------------------

        length = len(seq)

        # ------------------------------------------
        # Entropy
        # ------------------------------------------

        entropy = shannon_entropy(seq)

        # ------------------------------------------
        # Amino-acid composition
        # ------------------------------------------

        counts = Counter(seq)

        aa, aa_count = counts.most_common(1)[0]

        aa_percent = (
            aa_count / length
        ) * 100

        # ------------------------------------------
        # Mutation %
        # ------------------------------------------

        best_k = None

        best_mut = length

        for k in range(
            1,
            min(MAX_K, length) + 1
        ):

            mut = mutations_to_repeat(
                seq,
                k
            )

            if mut < best_mut:

                best_mut = mut

                best_k = k

        mutation_percent = (
            best_mut / length
        ) * 100

        rows.append({

            "Protein": protein,

            "Start": start,

            "End": end,

            "Length": length,

            "Entropy": round(
                entropy,
                4
            ),

            "Most_Common_AA": aa,

            "Most_Common_AA_Percent": round(
                aa_percent,
                2
            ),

            "Best_Model": best_k,

            "Mutation_Percent": round(
                mutation_percent,
                2
            )

        })

    df = pd.DataFrame(rows)

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
        f"Saved {outfile}"
    )

print("\nDone.")