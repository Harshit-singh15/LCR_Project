from Bio import SeqIO

# ==========================================================
# INPUT FILES
# ==========================================================

fasta_file = r"ecoli\ecoli_cleaned.fasta"

bed_file = r"ecoli\dataforFig8\ecoli_missing_residues.bed"

output_fasta = r"ecoli\dataforFig8\ecoli_missing_residues.fa"

log_file = r"ecoli\dataforFig8\extract_missing_sequences.log"

# ==========================================================
output_fasta.parent.mkdir(
    parents=True,
    exist_ok=True
)
log_file.parent.mkdir(
    parents=True,
    exist_ok=True
)

def get_accession(record):

    """
    Extract UniProt accession from FASTA header.

    Supports

    sp|P12345|NAME

    tr|Q9ABC1|NAME

    P12345
    """

    parts = record.id.split("|")

    if len(parts) >= 2:

        return parts[1]

    return record.id


print("Loading proteome...")

proteins = {}

for record in SeqIO.parse(fasta_file, "fasta"):

    proteins[
        get_accession(record)
    ] = str(record.seq)

print(f"Proteins loaded : {len(proteins):,}")

total_regions = 0
extracted = 0
missing_proteins = 0
invalid_coordinates = 0
empty_regions = 0

log = open(log_file, "w")

out = open(output_fasta, "w")

print("\nExtracting missing residue sequences...")

with open(bed_file) as f:

    for line in f:

        total_regions += 1

        acc, start, end = line.rstrip().split("\t")

        start = int(start)

        end = int(end)

        if acc not in proteins:

            missing_proteins += 1

            log.write(
                f"Protein not found : {acc}\n"
            )

            continue

        sequence = proteins[acc]

        if (
            start < 1
            or
            end > len(sequence)
            or
            start > end
        ):

            invalid_coordinates += 1

            log.write(
                f"Invalid coordinates : "
                f"{acc} "
                f"{start}-{end}\n"
            )

            continue

        region = sequence[start-1:end]

        if len(region) == 0:

            empty_regions += 1

            log.write(
                f"Empty region : "
                f"{acc} "
                f"{start}-{end}\n"
            )

            continue

        out.write(
            f">{acc}|{start}|{end}\n"
        )

        out.write(region + "\n")

        extracted += 1

out.close()

log.close()

print("\nDone.")

print("---------------------------")

print(f"BED regions           : {total_regions:,}")

print(f"Extracted             : {extracted:,}")

print(f"Missing proteins      : {missing_proteins:,}")

print(f"Invalid coordinates   : {invalid_coordinates:,}")

print(f"Empty regions         : {empty_regions:,}")

print("\nSaved:")

print(output_fasta)