# check_protein_ids.py

fasta_file = "outputs/mouse_proteome.fasta"
bed_file   = "bed_files\\use_reference_lcr.bed"

# -----------------------------
# Read FASTA IDs
# -----------------------------
fasta_ids = set()

with open(fasta_file) as f:
    for line in f:
        if line.startswith(">"):
            fasta_ids.add(line[1:].strip().split()[0])

# -----------------------------
# Read BED IDs
# -----------------------------
bed_ids = set()

with open(bed_file) as f:
    next(f)  # skip header

    for line in f:
        cols = line.rstrip().split("\t")

        if len(cols) < 3:
            continue

        bed_ids.add(cols[0])

# -----------------------------
# Compare
# -----------------------------
only_in_bed = bed_ids - fasta_ids
only_in_fasta = fasta_ids - bed_ids

print("=" * 60)
print("FASTA proteins       :", len(fasta_ids))
print("BED proteins         :", len(bed_ids))
print("Only in BED          :", len(only_in_bed))
print("Only in FASTA        :", len(only_in_fasta))
print("=" * 60)

if only_in_bed:
    print("\nFirst 20 IDs present in BED but not FASTA:")
    for x in sorted(list(only_in_bed))[:20]:
        print(x)

if only_in_fasta:
    print("\nFirst 20 IDs present in FASTA but not BED:")
    for x in sorted(list(only_in_fasta))[:20]:
        print(x)