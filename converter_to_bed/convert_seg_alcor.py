from pathlib import Path
import re

print("[INFO] Running FASTA to BED converter")

# 1. PASTE YOUR EXACT PATHS HERE
# Use absolute paths (e.g., r"C:\path\to\file.fa") or relative paths.
input_file = Path(r"mouse\lcrbytools\mouse_alcor_mode2_masked.fa")
output_file = Path(r"mouse\bed_mouse\alcor_mode2_masked_mouse.bed")

# Ensure the parent directory for the output file exists
output_file.parent.mkdir(parents=True, exist_ok=True)

# Validate that the input file exists before running
if not input_file.is_file():
    raise FileNotFoundError(f"Input file not found: {input_file}")

print(f"Processing: {input_file.name}")

# Process the single FASTA file
with output_file.open("w", encoding="utf-8") as out:
    current_id = None
    sequence = []

    with input_file.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                # Process the previous sequence block before moving to the new header
                if current_id is not None:
                    seq = "".join(sequence)
                    for match in re.finditer(r"[a-z]+", seq):
                        start = match.start() + 1
                        end = match.end()
                        out.write(f"{current_id}\t{start}\t{end}\n")

                # Parse the new header
                header = line[1:].strip()
                parts = header.split("|")

                if len(parts) >= 3:
                    current_id = parts[1]
                else:
                    current_id = header.split()[0]

                sequence = []
            else:
                sequence.append(line)

        # Process the final sequence block at the end of the file
        if current_id is not None:
            seq = "".join(sequence)
            for match in re.finditer(r"[a-z]+", seq):
                start = match.start() + 1
                end = match.end()
                out.write(f"{current_id}\t{start}\t{end}\n")

print(f"[SUCCESS] Output saved to: {output_file}")
print("Conversion complete.")