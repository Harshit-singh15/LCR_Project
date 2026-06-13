import re

# ===== EDIT THESE =====
input_fasta = r"outputs\flps_strict_masked.fa"
output_bed  = r"bed_files\flps_strict_masked_mouse.bed"
# ======================

with open(output_bed, "w") as out:

    out.write("Protein_ID\tStart\tEnd\n")

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

                current_id = line[1:].split()[0]
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

print("SEG conversion complete.")
print("Output:", output_bed)