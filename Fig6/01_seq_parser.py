from pathlib import Path
from Bio import SeqIO
import sys

input_fasta = Path(sys.argv[1])  # Proteome FASTA file
output_bed = Path(sys.argv[2])   # Output BED file

window_size = 20
step_size = 10

# Create the directory path directly
output_bed.parent.mkdir(
    parents=True,
    exist_ok=True
)

def main():
    if step_size <= 0 or window_size <= 0:
        raise ValueError(
            "Window and step sizes must be positive integers"
        )

    with open(output_bed, "w") as out:
        # Memory safe: streams the FASTA file one record at a time
        for record in SeqIO.parse(input_fasta, "fasta"):

            protein_id = (
                record.id.split("|")[1]
                if "|" in record.id
                else record.id
            )
            seq_len = len(record.seq)

            # FIX: Stop iterating when a full window can no longer fit
            for start in range(0, seq_len - window_size + 1, step_size):
                end = start + window_size

                # Convert to 1-based inclusive coordinates
                bed_start = start + 1
                bed_end = end

                # Memory safe: writes instantly to disk
                out.write(
                    f"{protein_id}\t{bed_start}\t{bed_end}\n"
                )

    print(f"Saved: {output_bed}")

if __name__ == "__main__":
    main()