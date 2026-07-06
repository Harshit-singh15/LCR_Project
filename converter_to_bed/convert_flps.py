from pathlib import Path
import sys


# Define your exact input file and desired output file paths
# (You can replace "target_file.out" and "target_file.bed" with your actual filenames)
input_file = Path(r"human\lcrbytools_human\human_flps_default.out")
output_file = Path(r"human\bed_human\human_flps_default.bed")

# Ensure the parent directory for the output file exists
output_file.parent.mkdir(parents=True, exist_ok=True)

# Validate that the input file exists before running
if not input_file.is_file():
    raise FileNotFoundError(f"Input file not found: {input_file}")

print(f"Processing: {input_file.name}")

written = 0
skipped = 0

# Process the single file
with output_file.open("w", encoding="utf-8") as out:
    with input_file.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            fields = line.split()

            # Skip malformed lines
            if len(fields) < 5:
                skipped += 1
                continue

            try:
                protein_field = fields[0]

                # UniProt accession parsing
                # tr|A0A067XG43|A0A067XG43_CAEEL -> A0A067XG43
                parts = protein_field.split("|")
                if len(parts) >= 3:
                    protein_id = parts[1]
                else:
                    protein_id = protein_field

                start = int(fields[3])
                end = int(fields[4])

                out.write(f"{protein_id}\t{start}\t{end}\n")
                written += 1

            except Exception:
                skipped += 1

print("\nFinished")
print(f"Output saved to: {output_file}")
print("Regions written:", written)
print("Skipped lines:", skipped)