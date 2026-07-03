"""
Figure 2 - Step 2
Extract FASTA sequences from consensus BED files.

Pure Python implementation of:

bedtools getfasta \
    -fi proteome.fasta \
    -bed consensus_x.bed \
    -name \
    -fo consensus_x.fa
"""

from pathlib import Path
import time


# ==========================================================
# Read FASTA once
# ==========================================================

def read_fasta(fasta_file):
    """
    Reads a protein FASTA file into a dictionary.

    Returns
    -------
    {
        protein_id : sequence
    }
    """

    sequences = {}

    header = None
    seq = []

    with open(fasta_file, "r") as f:

        for line in f:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if header is not None:
                    sequences[header] = "".join(seq)

                header = line[1:].split()[0]
                seq = []

            else:
                seq.append(line)

        if header is not None:
            sequences[header] = "".join(seq)

    return sequences


# ==========================================================
# Figure 2 - Step 2
# ==========================================================
from pathlib import Path
import time


# ==========================================================
# Read FASTA
# ==========================================================

def read_fasta(fasta_file):
    """
    Reads a protein FASTA file into a dictionary.

    Returns
    -------
    {
        protein_id : sequence
    }
    """

    fasta = {}

    header = None
    seq = []

    with open(fasta_file, "r") as f:

        for line in f:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if header is not None:
                    fasta[header] = "".join(seq)

                header = line[1:].split()[0]
                seq = []

            else:
                seq.append(line)

        if header is not None:
            fasta[header] = "".join(seq)

    return fasta


# ==========================================================
# Extract FASTA from BED
# ==========================================================

def extract_fasta_from_beds(
    fasta_file,
    bed_dir,
    output_dir,
    pattern="*.bed"
):
    """
    Extract FASTA sequences corresponding to BED regions.

    Parameters
    ----------
    fasta_file : str or Path
        Protein FASTA file.

    bed_dir : str or Path
        Folder containing BED files.

    output_dir : str or Path
        Folder where FASTA files will be written.

    pattern : str
        BED filename pattern.
        Examples:
            "*.bed"
            "consensus_*.bed"
    """

    overall_start = time.perf_counter()

    project_dir = Path.cwd()

    fasta_file = Path(fasta_file)
    bed_dir = Path(bed_dir)
    output_dir = Path(output_dir)

    if not fasta_file.is_absolute():
        fasta_file = project_dir / fasta_file

    if not bed_dir.is_absolute():
        bed_dir = project_dir / bed_dir

    if not output_dir.is_absolute():
        output_dir = project_dir / output_dir

    if not fasta_file.exists():
        raise FileNotFoundError(
            f"\nFASTA file not found:\n{fasta_file}"
        )

    if not bed_dir.exists():
        raise FileNotFoundError(
            f"\nBED folder not found:\n{bed_dir}"
        )

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    bed_files = sorted(
        bed_dir.glob(pattern)
    )

    if len(bed_files) == 0:
        raise RuntimeError(
            f"\nNo BED files found matching:\n{pattern}"
        )

    print("\n" + "=" * 60)
    print("Extract FASTA from BED")
    print("=" * 60)

    print("\nReading proteome...")

    sequences = read_fasta(fasta_file)

    print(f"Loaded {len(sequences):,} proteins.\n")

    total_sequences = 0

    processed = 0

    skipped = 0

    for bed_file in bed_files:

        if bed_file.stat().st_size == 0:

            print(f"Skipping empty file : {bed_file.name}")

            skipped += 1

            continue

        file_start = time.perf_counter()

        output_file = output_dir / f"{bed_file.stem}.fa"

        extracted = 0

        missing = 0

        with open(bed_file) as bed, open(output_file, "w") as out:

            for line in bed:

                line = line.strip()

                if not line:
                    continue

                cols = line.split("\t")

                if len(cols) < 3:
                    continue

                protein = cols[0]

                try:

                    start = int(cols[1])

                    end = int(cols[2])

                except ValueError:

                    continue

                if protein not in sequences:

                    missing += 1

                    continue

                seq = sequences[protein]

                start = max(0, start)

                end = min(end, len(seq))

                subseq = seq[start:end]

                if not subseq:
                    continue

                out.write(f">{protein}:{start}-{end}\n")
                out.write(subseq + "\n")

                extracted += 1

        processed += 1

        total_sequences += extracted

        elapsed = time.perf_counter() - file_start

        print(
            f"{bed_file.stem:<30}"
            f"{extracted:>7,} seqs"
            f"   {elapsed:>6.2f} sec"
        )

        if missing > 0:

            print(f"   Missing proteins : {missing}")

    total_time = time.perf_counter() - overall_start

    print("\n" + "=" * 60)
    print("Finished")
    print("=" * 60)

    print(f"BED files processed : {processed}")

    print(f"Skipped            : {skipped}")

    print(f"Sequences written  : {total_sequences:,}")

    print(f"Output folder      : {output_dir}")

    print(f"Total time         : {total_time:.2f} sec")

    print("=" * 60)

    return output_dir


extract_fasta_from_beds(
    fasta_file="ecoli\\ecoli_cleaned.fasta",
    bed_dir="ecoli\\dataforFig2\\02_consensus_beds",
    output_dir="ecoli\\dataforFig2\\03_consensus_fastas",
    pattern="consensus_*.bed"
)