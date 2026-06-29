import os

# ==================================
input_file = r"ecoli\lcrbytools\flps2_default_ecoli.out"
output_bed = r"ecoli\bed_ecoli\flps2_default_ecoli.bed"
# ==================================

os.makedirs(
    os.path.dirname(output_bed),
    exist_ok=True
)

written = 0
skipped = 0

with open(output_bed, "w") as out:

    with open(input_file, "r") as f:

        for line in f:

            line = line.strip()

            if not line:
                continue

            fields = line.split()

            # skip malformed lines
            if len(fields) < 5:
                skipped += 1
                continue

            try:

                protein_field = fields[0]

                # UniProt accession
                #
                # tr|A0A067XG43|A0A067XG43_CAEEL
                #
                # -> A0A067XG43

                parts = protein_field.split("|")

                if len(parts) >= 3:
                    protein_id = parts[1]
                else:
                    protein_id = protein_field

                start = int(fields[3])
                end   = int(fields[4])

                out.write(
                    f"{protein_id}\t{start}\t{end}\n"
                )

                written += 1

            except Exception:
                skipped += 1

print("Finished")
print("Regions written:", written)
print("Skipped lines:", skipped)
print("Output:", output_bed)