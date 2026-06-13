import pandas as pd

df = pd.read_csv("outputs/mouse_proteome_clustalw.tsv", sep="\t")

df[["seqid","start","end"]].rename(
    columns={
        "seqid":"Protein_ID",
        "start":"Start",
        "end":"End"
    }
).to_csv(
    "bed_files/treks_clustalw.bed",
    sep="\t",
    index=False
)