from pathlib import Path
import pandas as pd

# ==========================================================
# Directories
# ==========================================================

INPUT_DIR = Path(
    r"ecoli\dataforFig2\05_filtered"
)

OUTPUT_DIR = Path(
    r"ecoli\dataforFig2\06_peptide_counts"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

# ==========================================================
# Statistics
# ==========================================================

files_processed = 0
regions_processed = 0

# ==========================================================
# Process
# ==========================================================

for file in sorted(
    INPUT_DIR.glob("*_filtered.tsv")
):

    print(f"Processing {file.name}")

    df = pd.read_csv(
        file,
        sep="\t"
    )

    regions_processed += len(df)

    # --------------------------------------
    # Count peptide + type together
    # --------------------------------------

    counts = (
        df.groupby(
            ["Best-Peptide", "Best-Type"]
        )
        .size()
        .reset_index(name="Count")
    )

    counts = counts.sort_values(
        "Count",
        ascending=False
    )

    counts["Proportion"] = (
        counts["Count"] /
        counts["Count"].sum()
    )

    outfile = (
        OUTPUT_DIR /
        file.name.replace(
            "_filtered.tsv",
            "_peptide_counts.tsv"
        )
    )

    counts.to_csv(
        outfile,
        sep="\t",
        index=False
    )

    print(f"Saved: {outfile.name}")

    files_processed += 1

# ==========================================================
# Summary
# ==========================================================

print("\n====================================")

print(
    f"Files processed : {files_processed}"
)

print(
    f"Regions         : {regions_processed}"
)

print("====================================")