import re

input_fasta = r"outputs\flps_strict_masked.fa"
output_bed = r"bed_files\flps_strict_masked_mouse.bed"

with open(output_bed, "w") as out:

    seq_id = None
    seq = []

    with open(input_fasta) as f:
        for line in f:

            if line.startswith(">"):

                if seq_id is not None:

                    sequence = "".join(seq)

                    for match in re.finditer(r"[a-z]+|X+", sequence):
                        start = match.start() + 1
                        end = match.end()

                        out.write(
                             f"{seq_id}\t{start}\t{end}\n"
                        )

                seq_id = line[1:].strip().split()[0]
                seq = []

            else:
                seq.append(line.strip())

        # process last sequence
        if seq_id is not None:

            sequence = "".join(seq)

            for match in re.finditer(r"[a-z]+|X+", sequence):
                start = match.start() + 1
                end = match.end()

                out.write(
                    f"{seq_id}\t{start}\t{end}\n"
                )

print("BED file written:", output_bed)