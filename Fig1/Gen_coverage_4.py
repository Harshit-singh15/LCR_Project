from pathlib import Path
from Bio import SeqIO
import pandas as pd
import gc
import sys

# =====================================================
# CONFIG
# =====================================================

EXTRACT_DIR = Path(sys.argv[1])
PROTEIN_LENGTHS = Path(sys.argv[2])
OUTPUT_DIR = Path(sys.argv[3])
FASTA_FILE = Path(sys.argv[4])

Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

# =====================================================
# LOAD PROTEIN LENGTHS
# =====================================================

lengths_df = pd.read_csv(
    PROTEIN_LENGTHS,
    sep="\t",
    usecols=["Protein_ID"],
    dtype={"Protein_ID": "string"}
)

TOTAL_PROTEINS = len(lengths_df)

all_proteins = set(lengths_df["Protein_ID"])

# =====================================================
# LOAD FASTA LENGTH LOOKUP
# =====================================================

length_dict = {}

for record in SeqIO.parse(FASTA_FILE, "fasta"):

    full_id = record.description.split()[0]

    length = len(record.seq)

    length_dict[full_id] = length

    if "|" in full_id:

        parts = full_id.split("|")

        if len(parts) >= 3:

            accession = parts[1]
            entry_name = parts[2]

            length_dict[accession] = length
            length_dict[entry_name] = length

            if "_" in entry_name:

                short_name = entry_name.split("_")[0]

                length_dict[short_name] = length

print(f"Loaded {len(length_dict)} protein identifiers")

# =====================================================
# MERGE INTERVALS
# =====================================================

def merged_length(intervals):

    if not intervals:
        return 0

    intervals.sort()

    total = 0

    cur_start, cur_end = intervals[0]

    for start, end in intervals[1:]:

        if start <= cur_end + 1:

            if end > cur_end:
                cur_end = end

        else:

            total += cur_end - cur_start + 1

            cur_start, cur_end = start, end

    total += cur_end - cur_start + 1

    return total

# =====================================================
# PROCESS FILES
# =====================================================

for file in EXTRACT_DIR.glob("*_lcrs.tsv"):

    print(f"\nProcessing {file.name}")

    df = pd.read_csv(
        file,
        sep="\t",
        usecols=["Protein_ID", "Start", "End"],
        dtype={
            "Protein_ID": "string",
            "Start": "int32",
            "End": "int32"
        }
    )

    intervals = {}

    for row in df.itertuples(index=False):

        intervals.setdefault(
            row.Protein_ID,
            []
        ).append(
            (
                row.Start,
                row.End
            )
        )

    proteins_with_lcr = set(intervals.keys())

    missing_ids = 0

    bin_counts = {
        "0-20": 0,
        "20-40": 0,
        "40-60": 0,
        "60-80": 0,
        "80-100": 0
    }

    # ----------------------------------------------
    # Proteins with LCR
    # ----------------------------------------------

    for pid, ivals in intervals.items():

        if pid not in length_dict:

            print(f"Missing length: {pid}")

            missing_ids += 1

            continue

        coverage = (
            merged_length(ivals)
            /
            length_dict[pid]
        ) * 100

        if coverage < 20:

            bin_counts["0-20"] += 1

        elif coverage < 40:

            bin_counts["20-40"] += 1

        elif coverage < 60:

            bin_counts["40-60"] += 1

        elif coverage < 80:

            bin_counts["60-80"] += 1

        else:

            bin_counts["80-100"] += 1

    # ----------------------------------------------
    # Proteins without LCR
    # ----------------------------------------------

    bin_counts["0-20"] += (
        TOTAL_PROTEINS -
        len(proteins_with_lcr)
    )

    counts = [

        ["0-20", bin_counts["0-20"]],
        ["20-40", bin_counts["20-40"]],
        ["40-60", bin_counts["40-60"]],
        ["60-80", bin_counts["60-80"]],
        ["80-100", bin_counts["80-100"]]

    ]

    out_df = pd.DataFrame(
        counts,
        columns=[
            "Category",
            "Count"
        ]
    )

    total = out_df["Count"].sum()

    print(f"Proteins counted: {total}")

    print(f"Missing IDs: {missing_ids}")

    if total != TOTAL_PROTEINS:

        print(
            f"WARNING: expected {TOTAL_PROTEINS}, got {total}"
        )

    output_file = (
        OUTPUT_DIR
        /
        f"{file.stem.replace('_lcrs','')}_categorized.tsv"
    )

    out_df.to_csv(
        output_file,
        sep="\t",
        index=False
    )

    del df
    del intervals
    del proteins_with_lcr
    del out_df
    gc.collect()

print("\nCoverage calculation complete.")