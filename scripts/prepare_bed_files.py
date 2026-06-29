from pathlib import Path
import pandas as pd

# ==========================================================
# INPUT / OUTPUT 
# ==========================================================

INPUT_DIR = Path(r"ecoli\bed_ecoli")      # original BED files
OUTPUT_DIR = Path(r"ecoli\bed_bedtools_Ecoli")     # BEDTools compatible files

OUTPUT_DIR.mkdir(exist_ok=True)

summary = []

# ==========================================================
# Process each BED
# ==========================================================


for bed in sorted(INPUT_DIR.glob("*.bed")):

    print(f"Processing {bed.name}")

    df = pd.read_csv(
        bed,
        sep=r"\s+",
        header=None,
        usecols=[0,1,2],
        names=["Protein","Start","End"],
        engine="python"
    )

    original = len(df)

    # --------------------------
    # remove missing rows
    # --------------------------

    df = df.dropna()

    df["Start"] = df["Start"].astype(int)
    df["End"]   = df["End"].astype(int)

    # --------------------------
    # convert 1-based -> BED
    #
    # BED:
    # start = start-1
    # end unchanged
    # --------------------------

    df["Start"] = df["Start"] - 1

    # BED start cannot be negative

    df.loc[df["Start"] < 0, "Start"] = 0

    # --------------------------
    # remove invalid intervals
    # --------------------------

    invalid = (df["Start"] >= df["End"]).sum()

    df = df[df["Start"] < df["End"]]

    # --------------------------
    # sort
    # --------------------------

    df = df.sort_values(
        ["Protein","Start","End"]
    )

    # --------------------------
    # save
    # --------------------------

    outfile = OUTPUT_DIR / bed.name

    df.to_csv(
        outfile,
        sep="\t",
        header=False,
        index=False
    )

    summary.append([
        bed.name,
        original,
        invalid,
        len(df)
    ])

# ==========================================================
# Summary
# ==========================================================

summary = pd.DataFrame(
    summary,
    columns=[
        "BED_File",
        "Original_Intervals",
        "Removed_Invalid",
        "Final_Intervals"
    ]
)

summary.to_csv(
    OUTPUT_DIR/"bedtools_summary.tsv",
    sep="\t",
    index=False
)

print("\n====================================")
print(summary)
print("====================================")