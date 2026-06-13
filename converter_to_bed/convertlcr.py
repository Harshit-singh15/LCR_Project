import pandas as pd

df = pd.read_csv("C:/Users/91892/Desktop/Scripts_lcr/outputs/LCRFinder_protein_sorted.tsv", sep="\t")
df[["Protein_ID","Start","End"]].to_csv(
    "lcrfinder.bed",
    sep="\t",
    index=False
)