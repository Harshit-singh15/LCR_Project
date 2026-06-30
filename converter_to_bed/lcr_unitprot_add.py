import pandas as pd

# ==========================
fasta_file = r"outputs\mouse_proteome.fasta"
lcrfinder_file = r"bed_files\lcrfinder.bed"
output_file = r"bed_files\lcrfinder_uniprot.bed"
# ==========================

# Build mapping
id_map = {}

counter = 1

with open(fasta_file) as f:

    for line in f:

        if line.startswith(">"):

            uniprot_id = line[1:].split()[0]

            id_map[counter] = uniprot_id

            counter += 1

# Read LCRFinder BED
df = pd.read_csv(lcrfinder_file, sep="\t")

# Replace numeric IDs
df["Protein_ID"] = df["Protein_ID"].map(id_map)

# Save
df.to_csv(output_file, sep="\t", index=False)

print("Done.")
print("Mapped", len(df), "LCRs")