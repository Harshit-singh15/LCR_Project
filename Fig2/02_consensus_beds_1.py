"""
Generate consensus BED files from BEDTools multiinter output.

Input:
    Fig2/multiinter.tsv

Output:
    Fig2/02_consensus_beds/
        consensus_1.bed
        consensus_2.bed
        ...
        consensus_summary.tsv

"""

from pathlib import Path
import pandas as pd

# ==========================================================
# INPUT / OUTPUT
# ==========================================================

INPUT_FILE = Path(r"zebrafish\multiinter.tsv")

OUTPUT_DIR = Path(r"zebrafish\dataforFig2\02_consensus_beds")

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

# ==========================================================
# Load multiinter
# ==========================================================

print("Reading multiinter.tsv ...")

multi = pd.read_csv(
    INPUT_FILE,
    sep="\t",
    header=None
)

# ----------------------------------------------------------
# Validation
# ----------------------------------------------------------

if multi.shape[1] < 4:
    raise ValueError(
        "multiinter.tsv has fewer than 4 columns."
    )

multi.columns = (
    ["Protein", "Start", "End", "Consensus"] +
    [f"Extra_{i}" for i in range(5, multi.shape[1] + 1)]
)

max_consensus = int(
    multi["Consensus"].max()
)

print(
    f"Detected consensus levels: 1 to {max_consensus}"
)

# ==========================================================
# Generate BEDs
# ==========================================================

summary = []

total_written = 0

for n in range(1, max_consensus + 1):

    subset = (
        multi[
            multi["Consensus"] == n
        ]
        [["Protein", "Start", "End"]]
        .copy()
    )

    outfile = (
        OUTPUT_DIR /
        f"consensus_{n}.bed"
    )

    subset.to_csv(
        outfile,
        sep="\t",
        header=False,
        index=False
    )

    count = len(subset)

    total_written += count

    summary.append([
        n,
        count,
        outfile.name
    ])

    print(
        f"Consensus {n:2d}: "
        f"{count:8,d} regions"
    )

# ==========================================================
# Save summary
# ==========================================================

summary_df = pd.DataFrame(
    summary,
    columns=[
        "Consensus",
        "Number_of_Regions",
        "BED_File"
    ]
)

summary_file = (
    OUTPUT_DIR /
    r"consensus_summary.tsv"
)

summary_df.to_csv(
    summary_file,
    sep="\t",
    index=False
)

# ==========================================================
# Final validation
# ==========================================================

if total_written != len(multi):

    raise RuntimeError(
        "Validation failed!\n"
        f"Rows in multiinter.tsv : {len(multi)}\n"
        f"Rows written          : {total_written}"
    )

print("\n=======================================")
print("Consensus BED generation completed.")
print(f"Total regions : {total_written:,}")
print(f"Summary saved : {summary_file}")
print("=======================================")