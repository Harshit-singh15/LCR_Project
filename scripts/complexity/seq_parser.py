from Bio import SeqIO


input_fasta = r"outputs\mouse_proteome.fasta"

window_size = 20
step_size = 10

output_bed = r"reference_lcr\mouse_windows.bed"

# ======================


def main():

    if step_size <= 0 or window_size <= 0:
        raise ValueError(
            "Window and step sizes must be positive integers"
        )

    with open(output_bed, "w") as out:

        for record in SeqIO.parse(input_fasta, "fasta"):

            protein_id = record.id
            seq_len = len(record.seq)

            for start in range(0, seq_len, step_size):

                end = start + window_size

                if start >= seq_len:
                    break

                end = min(end, seq_len)

                # Convert to 1-based inclusive coordinates
                bed_start = start + 1
                bed_end = end

                out.write(
                    f"{protein_id}\t{bed_start}\t{bed_end}\n"
                )

    print(f"Saved: {output_bed}")


if __name__ == "__main__":
    main()