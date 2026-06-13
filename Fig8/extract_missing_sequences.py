from Bio import SeqIO

# ==========================
# INPUT FILES
# ==========================

fasta_file = "outputs\\mouse_proteome.fasta"
bed_file   = "Fig8\\mouse_missing_residues.bed"

output_fasta = "Fig8\\mouse_missing_residues.fa"

# ==========================
# LOAD PROTEOME
# ==========================

proteins = {}

for record in SeqIO.parse(fasta_file, "fasta"):

    acc = record.id.split("|")[1]

    proteins[acc] = str(record.seq)

# ==========================
# EXTRACT REGIONS
# ==========================

out = open(output_fasta, "w")

count = 0

with open(bed_file) as f:

    for line in f:

        acc, start, end = line.strip().split("\t")

        start = int(start)
        end = int(end)

        if acc not in proteins:
            continue

        seq = proteins[acc]

        region = seq[start-1:end]

        if len(region) == 0:
            continue

        out.write(
            f">{acc}|{start}|{end}\n"
        )

        out.write(region + "\n")

        count += 1

out.close()

print(f"Extracted {count} regions")