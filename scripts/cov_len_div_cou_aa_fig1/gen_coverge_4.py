from pathlib import Path
import pandas as pd
from Bio import SeqIO

# =====================================================
# CONFIG
# =====================================================

EXTRACT_DIR = r"zebrafish\dataforFig1\extracted_sequences"
PROTEIN_LENGTHS = r"zebrafish\dataforFig1\protein_lengths.tsv"
OUTPUT_DIR = r"zebrafish\dataforFig1\LCR_Coverage"
FASTA_FILE = r"zebrafish\zebrafish.fasta"

TOTAL_PROTEINS = len(pd.read_csv(PROTEIN_LENGTHS, sep="\t"))

Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

# =====================================================
# LOAD PROTEIN LENGTHS


lengths_df = pd.read_csv(
    PROTEIN_LENGTHS,
    sep="\t"
)

# canonical proteins
all_proteins = set(
    lengths_df["Protein_ID"]
)

# lookup dictionary
length_dict = {}

for record in SeqIO.parse(
    FASTA_FILE,
    "fasta"
):

    full_id = record.description.split()[0]

    length = len(record.seq)

    # full id
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
# CANONICAL PROTEIN SET
# =====================================================

all_proteins = set(lengths_df["Protein_ID"])

# =====================================================
# INTERVAL MERGING
# =====================================================

def merged_length(intervals):

    if len(intervals) == 0:
        return 0

    intervals = sorted(intervals)

    merged = [list(intervals[0])]

    for start, end in intervals[1:]:

        last_start, last_end = merged[-1]

        if start <= last_end + 1:

            merged[-1][1] = max(
                last_end,
                end
            )

        else:

            merged.append(
                [start, end]
            )

    total = 0

    for start, end in merged:

        total += end - start + 1

    return total

# =====================================================
# PROCESS FILES
# =====================================================

for file in Path(EXTRACT_DIR).glob("*_lcrs.tsv"):

    print(f"\nProcessing {file.name}")

    df = pd.read_csv(
        file,
        sep="\t"
    )

    # ----------------------------------------------
    # Build interval list per protein
    # ----------------------------------------------

    intervals = {}

    for _, row in df.iterrows():

        pid = str(row["Protein_ID"])

        intervals.setdefault(
            pid,
            []
        ).append(
            (
                int(row["Start"]),
                int(row["End"])
            )
        )

    coverage_values = []

    missing_ids = 0

    # ----------------------------------------------
    # Proteins with LCRs
    # ----------------------------------------------

    for pid, ivals in intervals.items():

        if pid not in length_dict:

            print(f"Missing length: {pid}")
            missing_ids += 1
            continue

        protein_length = length_dict[pid]

        covered_residues = merged_length(
            ivals
        )

        coverage = (
            covered_residues /
            protein_length
        ) * 100

        coverage_values.append(
            coverage
        )

    # ----------------------------------------------
    # Proteins without LCRs
    # ----------------------------------------------

    proteins_with_lcr = set(intervals.keys())

    proteins_without_lcr = (
        all_proteins -
        proteins_with_lcr
    )

    coverage_values.extend(
        [0.0] *
        len(proteins_without_lcr)
    )

    coverage_series = pd.Series(
        coverage_values
    )

    # ----------------------------------------------
    # Bin coverage
    # ----------------------------------------------

    counts = []

    counts.append([
        "0-20",
        int(
            (
                (coverage_series >= 0) &
                (coverage_series < 20)
            ).sum()
        )
    ])

    counts.append([
        "20-40",
        int(
            (
                (coverage_series >= 20) &
                (coverage_series < 40)
            ).sum()
        )
    ])

    counts.append([
        "40-60",
        int(
            (
                (coverage_series >= 40) &
                (coverage_series < 60)
            ).sum()
        )
    ])

    counts.append([
        "60-80",
        int(
            (
                (coverage_series >= 60) &
                (coverage_series < 80)
            ).sum()
        )
    ])

    counts.append([
        "80-100",
        int(
            (
                (coverage_series >= 80) &
                (coverage_series <= 100)
            ).sum()
        )
    ])

    out_df = pd.DataFrame(
        counts,
        columns=[
            "Category",
            "Count"
        ]
    )

    # ----------------------------------------------
    # Validation
    # ----------------------------------------------

    total = out_df["Count"].sum()

    print(
        f"Proteins counted: {total}"
    )

    print(
        f"Missing IDs: {missing_ids}"
    )

    if total != TOTAL_PROTEINS:

        print(
            f"WARNING: expected {TOTAL_PROTEINS}, got {total}"
        )

    # ----------------------------------------------
    # Save
    # ----------------------------------------------

    output_file = (
        Path(OUTPUT_DIR) /
        f"{file.stem.replace('_lcrs','')}_categorized.tsv"
    )

    out_df.to_csv(
        output_file,
        sep="\t",
        index=False
    )

print("\nCoverage calculation complete.")