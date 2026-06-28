import os
import pandas as pd
from Bio import SeqIO

# =====================================================
# EDIT PATHS
# =====================================================

proteome_fasta = r"fruitfly\fruitfly_cleaned.fasta"

reference_bed = r"fruitfly\dataforFig6\fruitfly_windows_real.bed"

tool_folder = r"fruitfly\bed_bedtools_Fruitfly"

output_file = r"fruitfly\dataforFig6\protein_confusion.tsv"

# =====================================================


print("Loading proteome lengths...")

protein_lengths = {}

for record in SeqIO.parse(proteome_fasta, "fasta"):
    protein_lengths[record.id] = len(record.seq)

print(f"Loaded {len(protein_lengths):,} proteins")


# =====================================================
# REFERENCE LCR POSITIONS
# =====================================================

print("Loading reference LCR annotations...")

ref_df = pd.read_csv(
    reference_bed,
    sep=r"\s+",
    engine="python"
)

ref_df = ref_df[
    ref_df["Classification"] == "LCR"
]

reference_positions = {}

for protein in protein_lengths:
    reference_positions[protein] = set()

for _, row in ref_df.iterrows():

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
# PROCESS TOOLS
# =====================================================

results = []

tool_files = sorted(
    [
        f for f in os.listdir(tool_folder)
        if f.endswith(".bed")
    ]
)

print(f"Found {len(tool_files)} tool files")


for tool_file in tool_files:

    print(f"\nProcessing {tool_file}")

    tool_path = os.path.join(
        tool_folder,
        tool_file
    )

    tool_df = pd.read_csv(
        tool_path,
        sep=r"\s+",
        header=None,
        names=[
            "Protein_ID",
            "Start",
            "End"
        ],
        engine="python"
    )

    tool_positions = {}

    for protein in protein_lengths:
        tool_positions[protein] = set()

    for _, row in tool_df.iterrows():

        protein = row["Protein_ID"]

        start = int(row["Start"])
        end = int(row["End"])

        if protein not in tool_positions:
            tool_positions[protein] = set()

        tool_positions[protein].update(
            range(start, end + 1)
        )

    for protein in protein_lengths:

        length = protein_lengths[protein]

        ref = reference_positions.get(
            protein,
            set()
        )

        pred = tool_positions.get(
            protein,
            set()
        )

        tp = len(
            ref & pred
        )

        fp = len(
            pred - ref
        )

        fn = len(
            ref - pred
        )

        tn = (
            length
            - tp
            - fp
            - fn
        )

        if tn < 0:
            tn = 0

        results.append(
            [
                protein,
                tool_file.replace(".bed", ""),
                tp,
                fp,
                fn,
                tn
            ]
        )


# =====================================================
# SAVE
# =====================================================

out_df = pd.DataFrame(
    results,
    columns=[
        "Protein_ID",
        "Tool",
        "TP",
        "FP",
        "FN",
        "TN"
    ]
)

out_df.to_csv(
    output_file,
    sep="\t",
    index=False
)

print("\nSaved:")
print(output_file)

print(
    f"\nRows: {len(out_df):,}"
)