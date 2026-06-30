from Bio import SeqIO
from collections import Counter
from math import log2
import pandas as pd
import glob
import os

os.makedirs("Fig2/04_metrics", exist_ok=True)

rows = []

for fasta in glob.glob("Fig2/03_consensus_fastas/*.fa"):

    consensus = int(
        fasta.split("_")[-1].replace(".fa", "")
    )

    for record in SeqIO.parse(fasta, "fasta"):

        seq = str(record.seq)

        if len(seq) == 0:
            continue

        counts = Counter(seq)

        H = 0

        for c in counts.values():

            p = c / len(seq)

            H -= p * log2(p)

        rows.append([
            consensus,
            record.id,
            len(seq),
            H
        ])

entropy_df = pd.DataFrame(
    rows,
    columns=[
        "Consensus",
        "Sequence",
        "Length",
        "Entropy"
    ]
)

entropy_df.to_csv(
    "Fig2/04_metrics/entropy.tsv",
    sep="\t",
    index=False
)

print(entropy_df.head())