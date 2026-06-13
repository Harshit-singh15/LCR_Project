import pandas as pd
import os
import glob
import re

# ========= EDIT =========
input_dir = "bed_files"
output_dir = "bed_clean"
# ========================

os.makedirs(output_dir, exist_ok=True)

for file in glob.glob(os.path.join(input_dir, "*.bed")):

    print(f"Processing {os.path.basename(file)}")

    df = pd.read_csv(file, sep="\t")

    # Keep only first 3 columns
    df = df.iloc[:, :3]

    df.columns = ["Protein_ID", "Start", "End"]

    # Convert IDs
    def clean_id(pid):

        pid = str(pid)

        # UniProt style:
        # tr|A0A9L6KDR9|A0A9L6KDR9_MOUSE
        if "|" in pid:
            parts = pid.split("|")

            if len(parts) >= 2:
                return parts[1]

        return pid

    df["Protein_ID"] = df["Protein_ID"].apply(clean_id)

    # Ensure numeric
    df["Start"] = pd.to_numeric(df["Start"])
    df["End"] = pd.to_numeric(df["End"])

    # Sort
    df = df.sort_values(
        ["Protein_ID", "Start", "End"]
    )

    outfile = os.path.join(
        output_dir,
        os.path.basename(file)
    )

    df.to_csv(
        outfile,
        sep="\t",
        index=False,
        header=False
    )

print("Done")