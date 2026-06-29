import os
import pandas as pd

# ==================================
input_tsv = r"ecoli\lcrbytools\treks_clustalw_ecoli.tsv"
output_bed = r"ecoli\bed_ecoli\treks_clustalw_ecoli.bed"
# ==================================

os.makedirs(
    os.path.dirname(output_bed),
    exist_ok=True
)

df = pd.read_csv(input_tsv, sep="\t")

# Extract UniProt accession
def get_accession(seqid):

    parts = str(seqid).split("|")

    if len(parts) >= 3:
        return parts[1]

    return seqid

df["Protein_ID"] = df["seqid"].apply(get_accession)

df[["Protein_ID", "start", "end"]].rename(
    columns={
        "start": "Start",
        "end": "End"
    }
).to_csv(
    output_bed,
    sep="\t",
    index=False,
    header=False
)

print("T-REKS BED file written:", output_bed)
print("Regions:", len(df))