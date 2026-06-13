import pandas as pd
import re

# ==========================
# INPUT FILES
# ==========================

lcr_tsv = r"outputs\5_01_mouse_UP000000589_10090_21990_lcrbylcr.tsv"          # Zenodo TSV
proteome_fasta = r"outputs\mouse_proteome.fasta"

# ==========================
# OUTPUT FILES
# ==========================

bed_out = "bed_filesuse_reference_lcr.bed"
ambiguous_out = "ambiguous_lcrs.tsv"
unmapped_out = "unmapped_lcrs.tsv"

# ==========================
# LOAD PROTEOME
# ==========================

proteins = {}

current_id = None
seq_parts = []

with open(proteome_fasta) as f:
    for line in f:

        if line.startswith(">"):

            if current_id is not None:
                proteins[current_id] = "".join(seq_parts)

            # >tr|A0A1D5RM95|A0A1D5RM95_MOUSE ...
            fields = line.split("|")

            if len(fields) >= 2:
                current_id = fields[1]
            else:
                current_id = line[1:].split()[0]

            seq_parts = []

        else:
            seq_parts.append(line.strip())

    if current_id is not None:
        proteins[current_id] = "".join(seq_parts)

print(f"Loaded {len(proteins)} proteins")

# ==========================
# LOAD LCR DATA
# ==========================

df = pd.read_csv(lcr_tsv, sep="\t", index_col=0)

# Expected columns:
# parent
# sequence
# length
# LCR
# Species

bed_records = []
ambiguous_records = []
unmapped_records = []

# ==========================
# MAP LCRs
# ==========================

for lcr_id, row in df.iterrows():

    parent = str(row["parent"]).strip()
    lcr_seq = str(row["sequence"]).strip()

    if parent not in proteins:

        unmapped_records.append(
            [lcr_id, parent, "parent_not_found"]
        )
        continue

    protein_seq = proteins[parent]

    matches = [
        m.start()
        for m in re.finditer(
            re.escape(lcr_seq),
            protein_seq
        )
    ]

    if len(matches) == 1:

        start = matches[0] + 1      # 1-based
        end = start + len(lcr_seq) - 1

        bed_records.append(
            [parent, start, end, lcr_id]
        )

    elif len(matches) > 1:

        ambiguous_records.append(
            [lcr_id, parent, len(matches)]
        )

    else:

        unmapped_records.append(
            [lcr_id, parent, "sequence_not_found"]
        )

# ==========================
# WRITE OUTPUTS
# ==========================

pd.DataFrame(
    bed_records,
    columns=["Protein_ID", "Start", "End", "LCR_ID"]
).to_csv(
    bed_out,
    sep="\t",
    index=False
)

pd.DataFrame(
    ambiguous_records,
    columns=["LCR_ID", "Parent", "Num_Matches"]
).to_csv(
    ambiguous_out,
    sep="\t",
    index=False
)

pd.DataFrame(
    unmapped_records,
    columns=["LCR_ID", "Parent", "Reason"]
).to_csv(
    unmapped_out,
    sep="\t",
    index=False
)

print("\nFinished.")
print("BED:", len(bed_records))
print("Ambiguous:", len(ambiguous_records))
print("Unmapped:", len(unmapped_records))