from Bio import SeqIO
from collections import Counter
import pandas as pd
from pathlib import Path
import glob
from pathlib import Path

INPUT_DIR = Path("Fig3/01_fastas")

print("INPUT_DIR =", INPUT_DIR.resolve())
print("Exists:", INPUT_DIR.exists())

files = list(INPUT_DIR.glob("*.fa"))

print("Number of FASTA files:", len(files))

for f in files[:5]:
    print(f)

INPUT_DIR = Path(
    "Fig3\\01_fastas"
)

OUTPUT_DIR = Path(
    "Fig3\\02_metrics"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

for fasta in INPUT_DIR.glob("*.fa"):

    method = fasta.stem

    rows = []

    for record in SeqIO.parse(
        fasta,
        "fasta"
    ):

        seq = str(record.seq)

        if len(seq) == 0:
            continue

        counts = Counter(seq)

        purity = (
            max(counts.values())
            / len(seq)
        )

        rows.append(
            [purity]
        )

    df = pd.DataFrame(
        rows,
        columns=["Purity"]
    )

    outfile = (
        OUTPUT_DIR /
        f"{method}_purity.tsv"
    )

    df.to_csv(
        outfile,
        sep="\t",
        index=False
    )

    print(
        f"Saved {outfile}"
    ) 