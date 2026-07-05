from pathlib import Path

import pandas as pd
import sys
# ======================================================
# EDIT THESE
# ======================================================

input_file = Path(sys.argv[1])  # Input file with metrics
output_file = Path(sys.argv[2])  # Output file for classified regions

output_file.parent.mkdir(
    parents=True,
    exist_ok=True
)
# ======================================================

PRIORITY = {
    "Unknown": -1,
    "HCR": 0,
    "CBR": 1,
    "LCR": 2
}


def merge_intervals(windows):

    boundaries = set()

    for s, e, _ in windows:
        boundaries.add(s)
        boundaries.add(e + 1)

    boundaries = sorted(boundaries)

    merged = []

    for i in range(len(boundaries) - 1):

        seg_start = boundaries[i]
        seg_end = boundaries[i + 1] - 1

        best = None
        best_priority = -1

        for w_start, w_end, cls in windows:

            if w_start <= seg_start and w_end >= seg_end:

                p = PRIORITY[cls]

                if p > best_priority:
                    best_priority = p
                    best = cls

        if best is None:
            continue

        if (
            merged
            and merged[-1][2] == best
            and merged[-1][1] + 1 == seg_start
        ):
            merged[-1][1] = seg_end
        else:
            merged.append(
                [seg_start, seg_end, best]
            )

    return merged


def main():

    print("Reading file...")

    df = pd.read_csv(
        input_file,
        sep="\t"
    )

    df = df.sort_values(
        [
            "Protein_ID",
            "Start",
            "End"
        ]
    )

    total_proteins = df["Protein_ID"].nunique()

    print(f"Proteins : {total_proteins}")

    with open(output_file, "w") as out:

        out.write(
            "Protein_ID\tStart_Position\tEnd_Position\tClassification\n"
        )

        current_protein = None
        windows = []

        processed = 0

        for row in df.itertuples(index=False):

            protein = row.Protein_ID

            if current_protein is None:
                current_protein = protein

            if protein != current_protein:

                merged = merge_intervals(windows)

                for start, end, cls in merged:
                    out.write(
                        f"{current_protein}\t{start}\t{end}\t{cls}\n"
                    )

                processed += 1

                if processed % 500 == 0 or processed == total_proteins:
                    print(f"{processed}/{total_proteins}")

                windows = []

                current_protein = protein

            windows.append(
                (
                    int(row.Start),
                    int(row.End),
                    row.Classification
                )
            )

        if windows:

            merged = merge_intervals(windows)

            for start, end, cls in merged:
                out.write(
                    f"{current_protein}\t{start}\t{end}\t{cls}\n"
                )

            processed += 1

            print(f"{processed}/{total_proteins}")

    print()
    print("Done.")
    print("Saved:", output_file)


if __name__ == "__main__":
    main()