import pandas as pd
import re
import os

# ==========================
# INPUT FILES
# ==========================

lcr_tsv = r"arabidopsis\5_05_arabidopsis_UP000006548_3702_27468_lcrbylcr.tsv"
proteome_fasta = r"arabidopsis\arabidopsis.fasta"

# ==========================
# OUTPUT FILES
# ==========================

bed_out = r"arabidopsis\bed_arabidopsis\arabidopsis_dotplot.bed"
unmapped_out = r"arabidopsis\unmapped_lcrs_Arabidopsis.tsv"

# ==========================
# CREATE OUTPUT DIRECTORIES
# ==========================

os.makedirs(os.path.dirname(bed_out), exist_ok=True)
os.makedirs(os.path.dirname(unmapped_out), exist_ok=True)

# ==========================
# LOAD PROTEOME
# ==========================

proteins = {}

current_header = None
seq_parts = []

with open(proteome_fasta) as f:

    for line in f:

        if line.startswith(">"):

            if current_header is not None:

                seq = "".join(seq_parts)

                fields = current_header.split("|")

                if len(fields) >= 3:

                    accession = fields[1]
                    entry_name = fields[2].split()[0]

                    proteins[accession] = seq
                    proteins[entry_name] = seq

                    # SRA17_CAEEL -> SRA17
                    short_name = entry_name.split("_")[0]
                    proteins[short_name] = seq

                else:

                    protein_id = current_header[1:].split()[0]
                    proteins[protein_id] = seq

            current_header = line.strip()
            seq_parts = []

        else:
            seq_parts.append(line.strip())

    # last sequence
    if current_header is not None:

        seq = "".join(seq_parts)

        fields = current_header.split("|")

        if len(fields) >= 3:

            accession = fields[1]
            entry_name = fields[2].split()[0]

            proteins[accession] = seq
            proteins[entry_name] = seq

            short_name = entry_name.split("_")[0]
            proteins[short_name] = seq

        else:

            protein_id = current_header[1:].split()[0]
            proteins[protein_id] = seq

print(f"Loaded {len(proteins)} identifiers")

# ==========================
# LOAD LCR DATA
# ==========================

df = pd.read_csv(lcr_tsv, sep="\t", index_col=0)

bed_records = []
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

    if len(matches) > 0:

        for pos in matches:

            start = pos + 1      # 1-based
            end = start + len(lcr_seq) - 1

            bed_records.append(
                [parent, start, end]
            )

    else:

        unmapped_records.append(
            [lcr_id, parent, "sequence_not_found"]
        )

# ==========================
# WRITE BED
# ==========================

pd.DataFrame(
    bed_records,
    columns=["Protein_ID", "Start", "End"]
).to_csv(
    bed_out,
    sep="\t",
    index=False,
    header=False
)

# ==========================
# WRITE UNMAPPED
# ==========================

pd.DataFrame(
    unmapped_records,
    columns=["LCR_ID", "Parent", "Reason"]
).to_csv(
    unmapped_out,
    sep="\t",
    index=False
)

print("\nFinished")
print("BED entries:", len(bed_records))
print("Unmapped:", len(unmapped_records))