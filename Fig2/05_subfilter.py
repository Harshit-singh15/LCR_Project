from pathlib import Path
import pandas as pd

# ==========================================================
# Directories
# ==========================================================

INPUT_DIR = Path(r"ecoli\dataforFig2\04_substring_motifs")

OUTPUT_DIR = Path(r"ecoli\dataforFig2\05_filtered")

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
# Representative motif selection
# ==========================================================

def select_representative_motif(row):

    global regions_with_score_gt1

    candidates = [
        (
            "Mono",
            row["Mono-peptide"],
            row["Mono-Coverage"]
        ),
        (
            "Di",
            row["Di-peptide"],
            row["Di-Coverage"]
        ),
        (
            "Tri",
            row["Tri-peptide"],
            row["Tri-Coverage"]
        )
    ]

    # Keep only scores <= 1
    valid = []

    for motif_type, peptide, score in candidates:

        if score <= 1:

            valid.append(
                (
                    motif_type,
                    peptide,
                    score
                )
            )

        else:

            regions_with_score_gt1 += 1

    # Should never happen, but be safe
    if len(valid) == 0:

        return pd.Series(
            [
                "NA",
                "NA",
                0
            ]
        )

    # Highest score wins
    # If tie:
    # Tri > Di > Mono
    priority = {
        "Mono": 1,
        "Di": 2,
        "Tri": 3
    }

    best = max(
        valid,
        key=lambda x: (
            x[2],
            priority[x[0]]
        )
    )

    return pd.Series(
        [
            best[0],
            best[1],
            best[2]
        ]
    )

# ==========================================================
# Process one file
# ==========================================================

for file in sorted(
    INPUT_DIR.glob("*_substrings.tsv")
):

    print(
        f"Processing {file.name}"
    )

    df = pd.read_csv(
        file,
        sep="\t"
    )

    regions_processed += len(df)

    df[
        [
            "Best-Type",
            "Best-Peptide",
            "Best-Coverage"
        ]
    ] = df.apply(
        select_representative_motif,
        axis=1
    )

    result = df[
        [
            "Protein",
            "Start",
            "End",
            "Best-Type",
            "Best-Peptide",
            "Best-Coverage"
        ]
    ]

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

    print(
        f"Saved: {outfile.name}"
    )

    files_processed += 1

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