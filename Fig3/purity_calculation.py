from pathlib import Path
from Bio import SeqIO
from collections import Counter
import pandas as pd
import re

# ==========================================
# Directories
# ==========================================

INPUT_DIR = Path(r"fruitfly\dataforFig3\01_fastas")
OUTPUT_DIR = Path(r"fruitfly\dataforFig3\02_metrics")

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

# ==========================================

header_pattern = re.compile(
    r"(.+):(\d+)-(\d+)"
)

for fasta in sorted(INPUT_DIR.glob("*.fa")):

    method = fasta.stem

    rows = []

    print(f"Processing {method}")

    for record in SeqIO.parse(fasta, "fasta"):

        seq = str(record.seq).upper()

        if len(seq) == 0:
            continue

        match = header_pattern.match(record.id)

        if match:

            protein, start, end = match.groups()

        else:

            protein = record.id
            start = ""
            end = ""

        counts = Counter(seq)

        aa, aa_count = counts.most_common(1)[0]

        purity = aa_count / len(seq)

        rows.append({

            "Protein": protein,

            "Start": start,

            "End": end,

            "Length": len(seq),

            "Most_Common_AA": aa,

            "AA_Count": aa_count,

            "Purity": purity

        })

    df = pd.DataFrame(rows)

    outfile = OUTPUT_DIR / f"{method}_purity.tsv"

    df.to_csv(
        outfile,
        sep="\t",
        index=False
    )

    print(f"Saved {outfile}")

print("\nFinished.")