import pandas as pd
from pathlib import Path

# -----------------------------
# Directories
# -----------------------------

INPUT_DIR = Path(
    "Fig2/04_metrics/substrings"
)

OUTPUT_DIR = Path(
    "Fig2/04_metrics/filtered"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

# -----------------------------
# Function
# -----------------------------

def filter_best_coverage(
    file_path
):

    df = pd.read_csv(
        file_path,
        sep="\t"
    )

    def select_best(row):

        candidates = {

            "Mono": (
                row["Mono-peptide"],
                row["Mono-Coverage"]
            ),

            "Di": (
                row["Di-peptide"],
                row["Di-Coverage"]
            ),

            "Tri": (
                row["Tri-peptide"],
                row["Tri-Coverage"]
            )
        }

        best_type, (
            best_peptide,
            best_coverage
        ) = max(

            candidates.items(),

            key=lambda x:
            (
                x[1][1]
                if x[1][1] <= 1
                else -1
            )
        )

        return pd.Series(
            [
                best_type,
                best_peptide,
                best_coverage
            ]
        )

    df[
        [
            "Best-Type",
            "Best-Peptide",
            "Best-Coverage"
        ]
    ] = df.apply(
        select_best,
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

    out_file = (
        OUTPUT_DIR /
        file_path.name.replace(
            "_substrings.tsv",
            "_filtered.tsv"
        )
    )

    result.to_csv(
        out_file,
        sep="\t",
        index=False
    )

    print(
        f"Saved: {out_file}"
    )

# -----------------------------
# Run
# -----------------------------

for file in INPUT_DIR.glob(
    "*_substrings.tsv"
):

    print(
        f"Processing {file.name}"
    )

    filter_best_coverage(
        file
    )

print("Done.")