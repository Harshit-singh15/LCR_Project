import csv
import os
import sys
from pathlib import Path
from Bio import SeqIO

# =====================================================
# INPUTS
# =====================================================

proteome_fasta = Path(sys.argv[1])
reference_bed = Path(sys.argv[2])
tool_folder = Path(sys.argv[3])
output_file = Path(sys.argv[4])

output_file.parent.mkdir(
    parents=True,
    exist_ok=True
)

# =====================================================
# LOAD PROTEIN LENGTHS
# =====================================================

print("Loading proteome lengths...")

protein_lengths = {
    record.id: len(record.seq)
    for record in SeqIO.parse(proteome_fasta, "fasta")
}

print(f"Loaded {len(protein_lengths):,} proteins")

# =====================================================
# LOAD REFERENCE LCRS
# =====================================================

print("Loading reference LCR annotations...")

reference_positions = {}

with open(reference_bed, newline="") as f:

    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        if row["Classification"] != "LCR":
            continue

        protein = row["Protein_ID"]

        start = int(row["Start_Position"])
        end = int(row["End_Position"])

        if protein not in reference_positions:
            reference_positions[protein] = set()

        reference_positions[protein].update(
            range(start, end + 1)
        )

print("Reference positions loaded")

# =====================================================
# TOOL FILES
# =====================================================

tool_files = sorted(
    f for f in os.listdir(tool_folder)
    if f.endswith(".bed")
)

print(f"Found {len(tool_files)} tool files")

# =====================================================
# OUTPUT
# =====================================================

with open(output_file, "w", newline="") as out:

    writer = csv.writer(
        out,
        delimiter="\t"
    )

    writer.writerow(
        [
            "Protein_ID",
            "Tool",
            "TP",
            "FP",
            "FN",
            "TN"
        ]
    )

    # =====================================================
    # PROCESS EACH TOOL
    # =====================================================

    for tool_file in tool_files:

        print(f"\nProcessing {tool_file}")

        tool_positions = {}

        tool_path = tool_folder / tool_file

        with open(tool_path) as f:

            reader = csv.reader(
                f,
                delimiter="\t"
            )

            for row in reader:

                if len(row) < 3:
                    continue

                protein = row[0]

                start = int(row[1])
                end = int(row[2])

                if protein not in tool_positions:
                    tool_positions[protein] = set()

                tool_positions[protein].update(
                    range(start, end + 1)
                )

        tool_name = tool_file.replace(".bed", "")

        for protein, length in protein_lengths.items():

            ref = reference_positions.get(
                protein,
                set()
            )

            pred = tool_positions.get(
                protein,
                set()
            )

            tp = len(ref & pred)
            fp = len(pred - ref)
            fn = len(ref - pred)

            tn = length - tp - fp - fn

            if tn < 0:
                tn = 0

            writer.writerow(
                [
                    protein,
                    tool_name,
                    tp,
                    fp,
                    fn,
                    tn
                ]
            )

print("\nSaved:")
print(output_file)