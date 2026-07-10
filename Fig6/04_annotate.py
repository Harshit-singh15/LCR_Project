from pathlib import Path
import pandas as pd
import sys
from collections import defaultdict

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
    # Sort the windows by Start and End position to ensure correct interval slicing
    windows.sort(key=lambda x: (x[0], x[1]))
    
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
    print("Reading file in chunks and grouping windows...")

    # Dictionary to collect windows: { protein_id: [(start, end, classification), ...] }
    protein_windows = defaultdict(list)
    
    # Process the file in chunks to keep memory footprint low
    chunk_size = 100000
    chunks = pd.read_csv(
        input_file,
        sep="\t",
        usecols=["Protein_ID", "Start", "End", "Classification"],
        chunksize=chunk_size
    )

    for chunk in chunks:
        for row in chunk.itertuples(index=False):
            protein_windows[row.Protein_ID].append(
                (int(row.Start), int(row.End), row.Classification)
            )

    total_proteins = len(protein_windows)
    print(f"Proteins collected: {total_proteins}")
    print("Processing and merging intervals...")

    with open(output_file, "w") as out:
        out.write(
            "Protein_ID\tStart_Position\tEnd_Position\tClassification\n"
        )

        processed = 0
        
        # Iterate through each protein's collected windows
        for protein, windows in protein_windows.items():
            # Pass the windows to be sorted and merged
            merged = merge_intervals(windows)

            for start, end, cls in merged:
                out.write(
                    f"{protein}\t{start}\t{end}\t{cls}\n"
                )

            processed += 1
            if processed % 500 == 0 or processed == total_proteins:
                print(f"{processed}/{total_proteins}")

    print()
    print("Done.")
    print("Saved:", output_file)


if __name__ == "__main__":
    main()