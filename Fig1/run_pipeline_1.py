from Bio import SeqIO
from pathlib import Path
import sys

# =====================================================
# CONFIG
# =====================================================

FASTA_FILE = sys.argv[1]
BED_FOLDER = sys.argv[2]
OUTPUT_DIR = sys.argv[3]

# =====================================================
# CREATE OUTPUT DIRECTORIES
# =====================================================

Path(OUTPUT_DIR).mkdir(exist_ok=True)

EXTRACT_DIR = Path(OUTPUT_DIR) / "extracted_sequences"
EXTRACT_DIR.mkdir(exist_ok=True)

# =====================================================
# LOAD FASTA
# =====================================================

print("\nLoading proteome...")

proteins = {}
protein_lengths = {}
canonical_lengths = []

for record in SeqIO.parse(FASTA_FILE, "fasta"):

    full_id = record.description.split()[0]
    sequence = str(record.seq)
    length = len(sequence)

    proteins[full_id] = sequence
    protein_lengths[full_id] = length

    if "|" in full_id:
        protein_id = full_id.split("|")[1]
    else:
        protein_id = full_id

    canonical_lengths.append((protein_id, length))

    # ----------------------------------------
    # create aliases
    # ----------------------------------------

    if "|" in full_id:

        parts = full_id.split("|")

        if len(parts) >= 3:

            accession = parts[1]
            entry_name = parts[2]

            proteins[accession] = sequence
            proteins[entry_name] = sequence

            protein_lengths[accession] = length
            protein_lengths[entry_name] = length

            if "_" in entry_name:

                short_name = entry_name.split("_")[0]

                proteins[short_name] = sequence
                protein_lengths[short_name] = length

print(f"Loaded {len(protein_lengths)} identifiers")

# =====================================================
# SAVE PROTEIN LENGTHS
# =====================================================

length_file = Path(OUTPUT_DIR) / "protein_lengths.tsv"

with open(length_file, "w") as out:

    out.write("Protein_ID\tProtein_Length\n")

    for protein_id, length in canonical_lengths:

        out.write(f"{protein_id}\t{length}\n")

print("Saved protein_lengths.tsv")

# =====================================================
# PROCESS BED FILES
# =====================================================

bed_files = sorted(Path(BED_FOLDER).glob("*.bed"))

print(f"\nFound {len(bed_files)} BED files")

for bed_file in bed_files:

    print(f"\nProcessing {bed_file.name}")

    missing_set = set()

    missing_ids = 0
    invalid_coords = 0

    bad_start = 0
    bad_end = 0
    bad_order = 0

    output_file = EXTRACT_DIR / f"{bed_file.stem}_lcrs.tsv"

    with open(output_file, "w") as out:

        out.write(
            "Protein_ID\tStart\tEnd\tLength\tProtein_Length\tSequence\n"
        )

        with open(bed_file) as infile:

            for line in infile:

                line = line.strip()

                if not line:
                    continue

                parts = line.split()

                if len(parts) < 3:
                    continue

                protein_id = parts[0]

                try:
                    start = int(float(parts[1]))
                    end = int(float(parts[2]))
                except:
                    invalid_coords += 1
                    continue

                if protein_id not in proteins:

                    missing_ids += 1
                    missing_set.add(protein_id)
                    continue

                protein_seq = proteins[protein_id]
                protein_len = len(protein_seq)

                if start < 1:

                    bad_start += 1
                    continue

                if end > protein_len:

                    bad_end += 1
                    continue

                if start > end:

                    bad_order += 1
                    continue

                lcr_seq = protein_seq[start - 1:end]

                out.write(
                    f"{protein_id}\t"
                    f"{start}\t"
                    f"{end}\t"
                    f"{len(lcr_seq)}\t"
                    f"{protein_len}\t"
                    f"{lcr_seq}\n"
                )

    print("Unique missing IDs:", len(missing_set))

    if missing_set:

        with open(
            EXTRACT_DIR / f"{bed_file.stem}_missing_ids.txt",
            "w"
        ) as miss:

            for pid in sorted(missing_set):

                miss.write(pid + "\n")

    extracted_count = sum(1 for _ in open(output_file)) - 1

    print(f"Extracted LCRs : {extracted_count}")
    print(f"Missing IDs    : {missing_ids}")
    print(f"Invalid coords : {invalid_coords}")
    print("start<1 :", bad_start)
    print("end>len :", bad_end)
    print("start>end :", bad_order)

print("\nPipeline step 1 complete.")