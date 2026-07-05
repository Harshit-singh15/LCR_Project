import os
from Bio import SeqIO
import csv

# ==========================================================
# INPUT FILES
# ==========================================================

proteome_fasta = r"ecoli\ecoli_cleaned.fasta"

sifts_file = r"uniprot_segments_observed.tsv"

output_file = r"ecoli\dataforFig8\ecoli_sifts_observed.tsv"

output_file.parent.mkdir(
    parents=True,
    exist_ok=True
)

def get_accession(record_id):
    """
    Extract UniProt accession from FASTA header.

    Supports:
        sp|P12345|NAME
        tr|Q9XYZ1|NAME
        P12345
    """

    parts = record_id.split("|")

    if len(parts) >= 2:
        return parts[1]

    return record_id


print("Loading  proteome...")

proteins = set()

for record in SeqIO.parse(proteome_fasta, "fasta"):

    proteins.add(
        get_accession(record.id)
    )

print(f"Proteins loaded : {len(proteins):,}")

print("\nFiltering observed residue file...")

total_rows = 0
kept_rows = 0
kept_proteins = set()

with open(
    sifts_file,
    "r",
    newline="",
    encoding="utf-8"
) as infile, open(
    output_file,
    "w",
    newline="",
    encoding="utf-8"
) as outfile:

    reader = csv.DictReader(
        infile,
        delimiter="\t"
    )

    writer = csv.DictWriter(
        outfile,
        fieldnames=reader.fieldnames,
        delimiter="\t"
    )

    writer.writeheader()

    for row in reader:

        total_rows += 1

        accession = row["SP_PRIMARY"]

        if accession in proteins:

            writer.writerow(row)

            kept_rows += 1

            kept_proteins.add(accession)

print("\nDone.")

print(f"Total SIFTS rows      : {total_rows:,}")
print(f"Rows retained         : {kept_rows:,}")
print(f"Proteins represented  : {len(kept_proteins):,}")

print("\nSaved:")
print(output_file)