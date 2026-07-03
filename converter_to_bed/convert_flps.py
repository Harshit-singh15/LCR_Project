from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pipeline import config

input_fol = Path(config.LCRBYTOOLS_DIR)
outputfolder = Path(config.BED_DIR)

if not input_fol.is_dir():
    raise FileNotFoundError(f"Input folder not found: {input_fol}")

outputfolder.mkdir(parents=True, exist_ok=True)

written = 0
skipped = 0
processed_files = 0

for input_file in sorted(input_fol.iterdir()):
    if not input_file.is_file() or input_file.suffix != ".out":
        continue

    input_path = input_file
    output_bed = outputfolder / f"{input_file.stem}.bed"

    processed_files += 1
    file_written = 0
    file_skipped = 0

    with output_bed.open("w", encoding="utf-8") as out:
        with input_path.open("r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()

                if not line:
                    continue

                fields = line.split()

                # skip malformed lines
                if len(fields) < 5:
                    file_skipped += 1
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
                    end = int(fields[4])

                    out.write(f"{protein_id}\t{start}\t{end}\n")
                    file_written += 1

                except Exception:
                    file_skipped += 1

    written += file_written
    skipped += file_skipped

    print(f"Processed: {input_file.name}")
    print("Regions written:", file_written)
    print("Skipped lines:", file_skipped)
    print("Output:", output_bed)

print("Finished")
print("Files processed:", processed_files)
print("Total regions written:", written)
print("Total skipped lines:", skipped)