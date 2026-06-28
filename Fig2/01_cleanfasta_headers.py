from pathlib import Path

# ==========================================
INPUT_FASTA = r"arabidopsis\arabidopsis.fasta"

OUTPUT_FASTA = r"arabidopsis\arabidopsis_cleaned.fasta"
# ==========================================

converted = 0

with open(INPUT_FASTA) as fin, \
     open(OUTPUT_FASTA, "w") as fout:

    for line in fin:

        if line.startswith(">"):

            header = line[1:].strip()

            # UniProt format
            # tr|A0A067XG43|A0A067XG43_CAEEL ...

            if "|" in header:

                parts = header.split("|")

                if len(parts) >= 2:

                    protein = parts[1]

                else:

                    protein = header.split()[0]

            else:

                protein = header.split()[0]

            fout.write(f">{protein}\n")

            converted += 1

        else:

            fout.write(line)

print(f"Proteins processed : {converted}")
print(f"Saved : {OUTPUT_FASTA}")