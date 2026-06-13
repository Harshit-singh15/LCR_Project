import os
import re
import glob
from collections import Counter
from pathlib import Path

# -----------------------------
# Directories
# -----------------------------

INPUT_DIR = Path("Fig2/03_consensus_fastas")

OUTPUT_DIR = Path(
    "Fig2/04_metrics/substrings"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

# -----------------------------
# Functions
# -----------------------------

def find_most_frequent_substring(
    sequence,
    length
):
    substr_counts = Counter(
        sequence[i:i+length]
        for i in range(
            len(sequence)-length+1
        )
    )

    most_common = (
        substr_counts.most_common(1)
    )

    if most_common:
        return most_common[0]

    return "-", 0


def process_fasta_file(
    file_path
):

    out_file = (
        OUTPUT_DIR /
        f"{file_path.stem}_substrings.tsv"
    )

    with open(file_path) as infile, \
         open(out_file, "w") as outfile:

        outfile.write(
            "Protein\tStart\tEnd\t"
            "Mono-peptide\tMono-Coverage\t"
            "Di-peptide\tDi-Coverage\t"
            "Tri-peptide\tTri-Coverage\n"
        )

        sequence = ""
        header = ""

        for line in infile:

            line = line.strip()

            if line.startswith(">"):

                # process previous seq

                if sequence and header:

                    process_record(
                        header,
                        sequence,
                        outfile
                    )

                header = line
                sequence = ""

            else:

                sequence += line

        # last record

        if sequence and header:

            process_record(
                header,
                sequence,
                outfile
            )

    print(
        f"Saved: {out_file}"
    )


def process_record(
    header,
    sequence,
    outfile
):

    match = re.match(
        r">(.+):(\d+)-(\d+)",
        header
    )

    if not match:
        return

    protein, start, end = (
        match.groups()
    )

    mono, mono_count = (
        find_most_frequent_substring(
            sequence,
            1
        )
    )

    di, di_count = (
        find_most_frequent_substring(
            sequence,
            2
        )
    )

    tri, tri_count = (
        find_most_frequent_substring(
            sequence,
            3
        )
    )

    seq_len = len(sequence)

    mono_cov = (
        len(mono) * mono_count
    ) / seq_len

    di_cov = (
        len(di) * di_count
    ) / seq_len

    tri_cov = (
        len(tri) * tri_count
    ) / seq_len

    outfile.write(
        f"{protein}\t"
        f"{start}\t"
        f"{end}\t"
        f"{mono}\t"
        f"{mono_cov:.6f}\t"
        f"{di}\t"
        f"{di_cov:.6f}\t"
        f"{tri}\t"
        f"{tri_cov:.6f}\n"
    )

# -----------------------------
# Run
# -----------------------------

for fasta in INPUT_DIR.glob("*.fa"):

    print(
        f"Processing {fasta.name}"
    )

    process_fasta_file(
        fasta
    )

print("Done.")