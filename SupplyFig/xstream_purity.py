import os
import math
from collections import Counter
from Bio import SeqIO

# ===================================
# CHANGE ONLY THESE PATHS
# ===================================

FASTA_FILE = r"outputs\mouse_proteome.fasta"

BED_FOLDER = r"SupplyFig\SupplFig9\xstream_bed"

OUTPUT_FOLDER = r"SupplyFig\SupplFig9\complexity_xstream"

# ===================================

def read_bed(filepath):

    regions = []

    with open(filepath) as f:

        next(f)  # skip header

        for line in f:

            line = line.rstrip()

            if not line:
                continue

            parts = line.split("\t")

            if len(parts) < 3:
                continue

            try:

                raw_id = parts[0]

                if "|" in raw_id:
                    protein = raw_id.split("|")[1]
                else:
                    protein = raw_id

                start = int(parts[-2])
                end = int(parts[-1])

                regions.append(
                    (
                        protein,
                        start,
                        end
                    )
                )

            except Exception as e:
                print("BAD LINE:", line[:100])

    return regions

def load_fasta(fasta_file):

    seqs = {}

    for rec in SeqIO.parse(fasta_file, "fasta"):

        if "|" in rec.id:

            toks = rec.id.split("|")

            if len(toks) >= 2:
                accession = toks[1]
            else:
                accession = rec.id

        else:
            accession = rec.id

        seqs[accession] = str(rec.seq)

    return seqs


def most_freq_aa(seq):

    counts = Counter(seq)

    aa, cnt = counts.most_common(1)[0]

    return aa, 100 * cnt / len(seq)







# ===================================
# MAIN
# ===================================

os.makedirs(
    OUTPUT_FOLDER,
    exist_ok=True
)

print("Loading FASTA...")

proteins = load_fasta(
    FASTA_FILE
)

print(
    "Proteins loaded:",
    len(proteins)
)

bed_files = sorted([
    x
    for x in os.listdir(BED_FOLDER)
    if x.endswith(".bed")
])

print(
    "BED files found:",
    len(bed_files)
)

for bed in bed_files:

    print("\nProcessing:", bed)

    bed_path = os.path.join(
        BED_FOLDER,
        bed
    )

    regions = read_bed(
        bed_path
    )

    print("Regions:", len(regions))

    number = ''.join(
        c for c in bed
        if c.isdigit()
    )

    out_file = os.path.join(
        OUTPUT_FOLDER,
        f"complexity_{number}.tsv"
    )

    with open(out_file, "w") as out:

        out.write(
            "Protein_ID\tStart\tEnd\tMost_Frequent_AA_Percent\n"
        )

        written = 0

        for protein, start, end in regions:

            if protein not in proteins:
                continue

            seq = proteins[protein]

            subseq = seq[start - 1:end]

            if len(subseq) == 0:
                continue

            aa, aa_pct = most_freq_aa(
                subseq
            )

            out.write(
                f"{protein}\t"
                f"{start}\t"
                f"{end}\t"
                f"{aa_pct:.2f}\n"
            )

            written += 1

    print(
        f"Saved: {out_file}"
    )

    print(
        f"Written rows: {written}"
    )