import re

input_fasta = r"outputs\mouse_proteome_alcor2_masked.fa"
output_bed  = r"bed_files\mouse_proteome_alcor2.bed"

with open(output_bed, "w") as out:

    out.write("Protein_ID\tStart\tEnd\n")

    current_id = None
    sequence = []

    with open(input_fasta) as f:
        for line in f:

            if line.startswith(">"):

                if current_id is not None:

                    seq = "".join(sequence)

                    for m in re.finditer(r"[a-z]+", seq):
                        start = m.start() + 1      # 1-based
                        end   = m.end()

                        out.write(
                            f"{current_id}\t{start}\t{end}\n"
                        )

                current_id = line.split()[0][1:]
                sequence = []

            else:
                sequence.append(line.strip())

        # last protein
        if current_id is not None:

            seq = "".join(sequence)

            for m in re.finditer(r"[a-z]+", seq):
                start = m.start() + 1
                end   = m.end()

                out.write(
                    f"{current_id}\t{start}\t{end}\n"
                )

print("BED file written.")