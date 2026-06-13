from Bio import SeqIO
from collections import Counter
import pandas as pd
import glob

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

        purity = max(counts.values()) / len(seq)

        rows.append([
            consensus,
            record.id,
            len(seq),
            purity
        ])

purity_df = pd.DataFrame(
    rows,
    columns=[
        "Consensus",
        "Sequence",
        "Length",
        "Purity"
    ]
)

purity_df.to_csv(
    "Fig2/04_metrics/purity.tsv",
    sep="\t",
    index=False
)

print(purity_df.head())