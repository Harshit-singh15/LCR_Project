from pathlib import Path
import pandas as pd

rows = []

EXTRACT_DIR = Path("outputs_prerequisite/extracted_sequences")
ENTROPY_DIR = Path("outputs_prerequisite/Diversity")

for file in EXTRACT_DIR.glob("*_lcrs.tsv"):

    tool = file.stem.replace("_lcrs", "")

    df = pd.read_csv(file, sep="\t")

    entropy_file = ENTROPY_DIR / f"{tool}_SNS"

    mean_entropy = None

    if entropy_file.exists():

        ent = pd.read_csv(entropy_file, sep="\t")

        mean_entropy = ent["Shannon_Entropy"].mean()

    rows.append({

        "Tool": tool,

        "LCRs":
            len(df),

        "Proteins_with_LCR":
            df["Protein_ID"].nunique(),

        "Total_LCR_Residues":
            df["Length"].sum(),

        "Mean_LCR_Length":
            round(df["Length"].mean(), 2),

        "Mean_Entropy":
            round(mean_entropy, 3)
            if mean_entropy is not None
            else None
    })

summary = pd.DataFrame(rows)

summary = summary.sort_values(
    by="LCRs",
    ascending=False
)

print(summary)

summary.to_csv(
    "outputs/qc_summary.tsv",
    sep="\t",
    index=False
)

print("\nSaved: outputs/qc_summary.tsv")