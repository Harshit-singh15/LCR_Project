from pathlib import Path
from collections import Counter

# ==========================================================
# Directories
# ==========================================================

INPUT_DIR = Path(r"zebrafish\dataforFig2\03_consensus_fastas")

OUTPUT_DIR = Path(r"zebrafish\dataforFig2\04_substring_motifs")

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

# ==========================================================
# Statistics
# ==========================================================

total_sequences = 0
processed = 0
skipped = 0

# ==========================================================
# Functions
# ==========================================================

def find_most_frequent_substring(sequence, k):
    """
    Returns the most frequent substring of length k
    and its occurrence count.
    """

    if len(sequence) < k:
        return "-", 0

    counts = Counter(
        sequence[i:i+k]
        for i in range(len(sequence)-k+1)
    )

    motif, count = counts.most_common(1)[0]

    return motif, count


def parse_header(header):
    """
    Accepts either

    >Protein:0-20

    or

    >Protein::0-20
    """

    header = header.strip()[1:]

    header = header.replace("::", ":")

    try:

        protein, coords = header.rsplit(":", 1)

        start, end = coords.split("-")

        return protein, int(start), int(end)

    except Exception:

        return None


def process_record(header, sequence, outfile):

    global processed, skipped

    info = parse_header(header)

    if info is None:

        print(f"Warning: malformed header: {header}")

        skipped += 1

        return

    protein, start, end = info

    if len(sequence) == 0:

        print(f"Warning: empty sequence: {protein}")

        skipped += 1

        return

    mono, mono_count = find_most_frequent_substring(sequence,1)
    di, di_count = find_most_frequent_substring(sequence,2)
    tri, tri_count = find_most_frequent_substring(sequence,3)

    seq_len = len(sequence)

    mono_cov = (len(mono) * mono_count) / seq_len
    di_cov   = (len(di)   * di_count)   / seq_len
    tri_cov  = (len(tri)  * tri_count)  / seq_len

    outfile.write(
        f"{protein}\t"
        f"{start}\t"
        f"{end}\t"
        f"{mono}\t{mono_cov:.6f}\t"
        f"{di}\t{di_cov:.6f}\t"
        f"{tri}\t{tri_cov:.6f}\n"
    )

    processed += 1


def process_fasta(fasta_file):

    global total_sequences

    output = OUTPUT_DIR / (
        fasta_file.stem + "_substrings.tsv"
    )

    with open(fasta_file) as fin, \
         open(output,"w") as fout:

        fout.write(
            "Protein\tStart\tEnd\t"
            "Mono-peptide\tMono-Coverage\t"
            "Di-peptide\tDi-Coverage\t"
            "Tri-peptide\tTri-Coverage\n"
        )

        header = None
        sequence = []

        for line in fin:

            line = line.strip()

            if line.startswith(">"):

                if header is not None:

                    total_sequences += 1

                    process_record(
                        header,
                        "".join(sequence),
                        fout
                    )

                header = line
                sequence = []

            else:

                sequence.append(line)

        if header is not None:

            total_sequences += 1

            process_record(
                header,
                "".join(sequence),
                fout)

    print(f"Saved : {output.name}")

# ==========================================================
# Run
# ==========================================================

for fasta in sorted(INPUT_DIR.glob("*.fa")):

    print(f"Processing {fasta.name}")

    process_fasta(fasta)

print("\n====================================")
print("Finished.")
print(f"Sequences read      : {total_sequences}")
print(f"Successfully parsed : {processed}")
print(f"Skipped             : {skipped}")
print("====================================")