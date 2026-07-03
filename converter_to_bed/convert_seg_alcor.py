from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pipeline import config

input_fol = Path(config.LCRBYTOOLS_DIR)
outputfolder = Path(config.BED_DIR)

if not input_fol.is_dir():
    raise FileNotFoundError(f"Input folder not found: {input_fol}")

outputfolder.mkdir(parents=True, exist_ok=True)

for input_fasta in sorted(input_fol.iterdir()):
    if not input_fasta.is_file() or input_fasta.suffix != ".fa":
        continue

    input_path = input_fasta
    output_bed = outputfolder / f"{input_fasta.stem}.bed"

    with output_bed.open("w", encoding="utf-8") as out:
        current_id = None
        sequence = []

        with input_path.open("r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()

                if not line:
                    continue

                if line.startswith(">"):
                    if current_id is not None:
                        seq = "".join(sequence)

                        for match in re.finditer(r"[a-z]+", seq):
                            start = match.start() + 1
                            end = match.end()
                            out.write(f"{current_id}\t{start}\t{end}\n")

                    header = line[1:].strip()
                    parts = header.split("|")

                    if len(parts) >= 3:
                        current_id = parts[1]
                    else:
                        current_id = header.split()[0]

                    sequence = []
                else:
                    sequence.append(line)

            if current_id is not None:
                seq = "".join(sequence)

                for match in re.finditer(r"[a-z]+", seq):
                    start = match.start() + 1
                    end = match.end()
                    out.write(f"{current_id}\t{start}\t{end}\n")

    print(f"Processed: {input_fasta.name}")
    print("Output:", output_bed)

print("conversion complete.")