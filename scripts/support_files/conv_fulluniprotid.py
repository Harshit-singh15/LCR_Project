from Bio import SeqIO

# ===== EDIT THESE =====

proteome_fasta = r"outputs\mouse_proteome.fasta"

input_bed = r"bed_files\Dotplot.bed"

output_bed = r"Dotplot_fixed.bed"

# ======================


print("Building UniProt ID mapping...")

id_map = {}

for record in SeqIO.parse(proteome_fasta, "fasta"):

    header = record.id

    parts = header.split("|")

    if len(parts) >= 2:

        accession = parts[1]

        id_map[accession] = header

print(f"Loaded {len(id_map):,} mappings")


mapped = 0
unmapped = 0

with open(input_bed) as fin, open(output_bed, "w") as fout:

    for line in fin:

        parts = line.strip().split()

        if len(parts) < 3:
            continue

        short_id = parts[0]
        start = parts[1]
        end = parts[2]

        if short_id in id_map:

            full_id = id_map[short_id]

            fout.write(
                f"{full_id}\t{start}\t{end}\n"
            )

            mapped += 1

        else:

            unmapped += 1

print(f"Mapped:   {mapped:,}")
print(f"Unmapped: {unmapped:,}")

print(f"\nSaved: {output_bed}")