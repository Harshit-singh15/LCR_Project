from pathlib import Path
from collections import Counter
from math import log2
from Bio import SeqIO
import pandas as pd

# ==========================================================
# Directories
# ==========================================================

INPUT_DIR = Path(
    r"ecoli\dataforFig2\03_consensus_fastas"
)

OUTPUT_DIR = Path(
    r"ecoli\dataforFig2\entropy"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

# ==========================================================
# Statistics
# ==========================================================

rows = []

files_processed = 0
sequences_processed = 0
empty_sequences = 0

# ==========================================================
# Process
# ==========================================================

for fasta in sorted(INPUT_DIR.glob("*.fa")):

    consensus = int(
        fasta.stem.split("_")[-1]
    )

    print(
        f"Processing Consensus {consensus}"
    )

    for record in SeqIO.parse(fasta, "fasta"):

        sequence = str(record.seq)

        if len(sequence) == 0:

            empty_sequences += 1
            continue

        counts = Counter(sequence)

        entropy = 0

        for c in counts.values():

            p = c / len(sequence)

            entropy -= p * log2(p)

        rows.append(
            [
                consensus,
                record.id,
                len(sequence),
                entropy
            ]
        )

        sequences_processed += 1

    files_processed += 1

# ==========================================================
# Save
# ==========================================================

entropy_df = pd.DataFrame(
    rows,
    columns=[
        "Consensus",
        "Sequence",
        "Length",
        "Entropy"
    ]
)

outfile = OUTPUT_DIR / "entropy.tsv"

entropy_df.to_csv(
    outfile,
    sep="\t",
    index=False
)

# ==========================================================
# Summary
# ==========================================================

print("\n===================================")

print(
    f"Files processed      : {files_processed}"
)

print(
    f"Sequences processed  : {sequences_processed}"
)

print(
    f"Empty sequences      : {empty_sequences}"
)

print(
    f"Saved : {outfile}"
)

print("===================================")