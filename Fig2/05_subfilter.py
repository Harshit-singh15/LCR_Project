from pathlib import Path
import pandas as pd
import gc
import sys

# ==========================================================
# Directories
# ==========================================================

INPUT_DIR = Path(sys.argv[1])     # substring motifs
OUTPUT_DIR = Path(sys.argv[2])    # filtered motifs

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

# ==========================================================
# Statistics
# ==========================================================

files_processed = 0
regions_processed = 0
regions_with_score_gt1 = 0

# ==========================================================
# Priority (used for tie-breaking)
# ==========================================================

priority = {
    "Mono": 1,
    "Di": 2,
    "Tri": 3
}

# ==========================================================
# Process one file
# ==========================================================

for file in sorted(
    INPUT_DIR.glob("*_substrings.tsv")
):

    print(f"Processing {file.name}")

    df = pd.read_csv(
        file,
        sep="\t"
    )

    regions_processed += len(df)

    # ------------------------------------------------------
    # Rename columns so itertuples() creates valid attributes
    # ------------------------------------------------------

    df.rename(
        columns={
            "Mono-peptide": "Mono_peptide",
            "Mono-Coverage": "Mono_Coverage",
            "Di-peptide": "Di_Peptide",
            "Di-Coverage": "Di_Coverage",
            "Tri-peptide": "Tri_Peptide",
            "Tri-Coverage": "Tri_Coverage"
        },
        inplace=True
    )

    best_types = []
    best_peptides = []
    best_coverages = []

    for row in df.itertuples(index=False):

        candidates = []

        # --------------------------
        # Mono
        # --------------------------

        if row.Mono_Coverage <= 1:

            candidates.append(
                (
                    "Mono",
                    row.Mono_peptide,
                    row.Mono_Coverage
                )
            )

        else:

            regions_with_score_gt1 += 1

        # --------------------------
        # Di
        # --------------------------

        if row.Di_Coverage <= 1:

            candidates.append(
                (
                    "Di",
                    row.Di_Peptide,
                    row.Di_Coverage
                )
            )

        else:

            regions_with_score_gt1 += 1

        # --------------------------
        # Tri
        # --------------------------

        if row.Tri_Coverage <= 1:

            candidates.append(
                (
                    "Tri",
                    row.Tri_Peptide,
                    row.Tri_Coverage
                )
            )

        else:

            regions_with_score_gt1 += 1

        # --------------------------
        # No valid motif
        # --------------------------

        if not candidates:

            best_types.append("NA")
            best_peptides.append("NA")
            best_coverages.append(0)

            continue

        # --------------------------
        # Highest score wins
        # Tie: Tri > Di > Mono
        # --------------------------

        best = max(
            candidates,
            key=lambda x: (
                x[2],
                priority[x[0]]
            )
        )

        best_types.append(best[0])
        best_peptides.append(best[1])
        best_coverages.append(best[2])

    # ------------------------------------------------------
    # Output
    # ------------------------------------------------------

    result = pd.DataFrame({

        "Protein": df["Protein"],
        "Start": df["Start"],
        "End": df["End"],
        "Best-Type": best_types,
        "Best-Peptide": best_peptides,
        "Best-Coverage": best_coverages

    })

    outfile = (
        OUTPUT_DIR /
        file.name.replace(
            "_substrings.tsv",
            "_filtered.tsv"
        )
    )

    result.to_csv(
        outfile,
        sep="\t",
        index=False
    )

    print(f"Saved: {outfile.name}")

    files_processed += 1

    del df
    del result
    del best_types
    del best_peptides
    del best_coverages
    gc.collect()

# ==========================================================
# Summary
# ==========================================================

print("\n========================================")

print(
    f"Files processed           : {files_processed}"
)

print(
    f"Regions processed         : {regions_processed}"
)

print(
    f"Scores > 1 ignored        : {regions_with_score_gt1}"
)

print("========================================")