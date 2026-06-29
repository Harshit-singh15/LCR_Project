import os
import re

# ===== EDIT THESE =====
input_fasta = r"ecoli\lcrbytools\seg_strict_ecoli.fa"
output_bed  = r"ecoli\bed_ecoli\seg_strict_ecoli.bed"
# ======================

# Create output directory if needed
os.makedirs(
    os.path.dirname(output_bed),
    exist_ok=True
)

with open(output_bed, "w") as out:

    current_id = None
    sequence = []

    with open(input_fasta, "r") as f:

        for line in f:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                # Process previous protein
                if current_id is not None:

                    seq = "".join(sequence)

                    for match in re.finditer(r"[a-z]+", seq):

                        start = match.start() + 1
                        end = match.end()

                        out.write(
                            f"{current_id}\t{start}\t{end}\n"
                        )

                # Extract UniProt accession
                header = line[1:].strip()
                parts = header.split("|")

                if len(parts) >= 3:
                    current_id = parts[1]
                else:
                    current_id = header.split()[0]

                sequence = []

            else:
                sequence.append(line)

        # Process final protein
        if current_id is not None:

            seq = "".join(sequence)

            for match in re.finditer(r"[a-z]+", seq):

                start = match.start() + 1
                end = match.end()

                out.write(
                    f"{current_id}\t{start}\t{end}\n"
                )

print("conversion complete.")
print("Output:", output_bed)