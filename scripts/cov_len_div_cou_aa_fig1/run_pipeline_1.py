from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from Bio import SeqIO
import pandas as pd
from pipeline import config

FASTA_FILE = config.FASTA_FILE
BED_FOLDER = config.BED_DIR
OUTPUT_DIR = config.FIG1_DATA_DIR

print("[INFO] Running Run_pipeline_1")

if not FASTA_FILE.exists():
    raise FileNotFoundError(f"FASTA file not found: {FASTA_FILE}")

if not BED_FOLDER.exists():
    raise FileNotFoundError(f"BED directory not found: {BED_FOLDER}")

# =====================================================
# CREATE OUTPUT DIRECTORIES
# =====================================================

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

EXTRACT_DIR = config.EXTRACTED_SEQUENCES_DIR
EXTRACT_DIR.mkdir(parents=True, exist_ok=True)

# =====================================================
# LOAD FASTA
# =====================================================

print("[INFO] Loading FASTA")
print("\nLoading proteome...")

proteins = {}
protein_lengths = {}

for record in SeqIO.parse(FASTA_FILE, "fasta-pearson"):

    full_id = record.description.split()[0]

    sequence = str(record.seq)

    length = len(sequence)

    proteins[full_id] = sequence
    protein_lengths[full_id] = length

    # ----------------------------------------
    # create aliases
    # ----------------------------------------

    if "|" in full_id:

        parts = full_id.split("|")

        if len(parts) >= 3:

            accession = parts[1]

            entry_name = parts[2]

            proteins[entry_name] = sequence

            if "_" in entry_name:
                short_name = entry_name.split("_")[0]
                proteins[short_name] = sequence

            proteins[accession] = sequence
            proteins[entry_name] = sequence

            protein_lengths[accession] = length
            protein_lengths[entry_name] = length

print(f"Loaded {len(protein_lengths)} identifiers")

# =====================================================
# SAVE PROTEIN LENGTHS
# =====================================================

canonical_lengths = []

for record in SeqIO.parse(FASTA_FILE, "fasta-pearson"):

    full_id = record.description.split()[0]

    if "|" in full_id:
        protein_id = full_id.split("|")[1]
    else:
        protein_id = full_id

    canonical_lengths.append(
        [protein_id, len(record.seq)]
    )

length_df = pd.DataFrame(
    canonical_lengths,
    columns=[
        "Protein_ID",
        "Protein_Length"
    ]
)

length_df.to_csv(
    config.PROTEIN_LENGTHS_FILE,
    sep="\t",
    index=False
)

print("Saved protein_lengths.tsv")

# =====================================================
# PROCESS BED FILES
# =====================================================

bed_files = sorted(BED_FOLDER.glob(config.BED_GLOB_PATTERN))

print("[INFO] Reading BED files")
print(f"[INFO] Found {len(bed_files)} BED files")

if not bed_files:
    raise FileNotFoundError("No BED files found in the configured BED directory")

for bed_file in bed_files:
    missing_list = []
    print(f"[INFO] Processing {bed_file.name}")

    try:

        df = pd.read_csv(
            bed_file,
            sep=r"\s+",
            header=None,
            engine="python"
        )

    except Exception as e:

        print(f"Could not read {bed_file}")
        print(e)
        continue

    # ----------------------------------------
    # keep first 3 columns only
    # ----------------------------------------

    df = df.iloc[:, :3]

    df.columns = [
        "Protein_ID",
        "Start",
        "End"
    ]

    extracted = []

    missing_ids = 0
    invalid_coords = 0


    for _, row in df.iterrows():

        protein_id = str(row["Protein_ID"])

        if pd.isna(row["Start"]) or pd.isna(row["End"]):
            invalid_coords += 1
            continue

        start = int(float(row["Start"]))
        end = int(float(row["End"]))


        if protein_id not in proteins:

            missing_ids += 1
            missing_list.append(protein_id)
            continue

        protein_seq = proteins[protein_id]

        protein_len = len(protein_seq)

        # ----------------------------------------
        # coordinate validation
        # ----------------------------------------
        bad_start = 0
        bad_end = 0
        bad_order = 0

        if start < 1:

            bad_start += 1
            continue

        if end > protein_len:

            bad_end += 1
            continue

        if start > end:

            bad_order += 1
            continue

        # ----------------------------------------
        # 0-based 
        # ----------------------------------------

        lcr_seq = protein_seq[start -1 :end]

        extracted.append(
            [
                protein_id,
                start,
                end,
                len(lcr_seq),
                protein_len,
                lcr_seq
            ]
        )

    out_df = pd.DataFrame(
        extracted,
        columns=[
            "Protein_ID",
            "Start",
            "End",
            "Length",
            "Protein_Length",
            "Sequence"
        ]
    )

    print("[INFO] Writing extracted sequences")

    output_file = (
        EXTRACT_DIR /
        f"{bed_file.stem}_lcrs.tsv"
    )

    out_df.to_csv(
        output_file,
        sep="\t",
        index=False
    )
    print("Unique missing IDs:", len(set(missing_list)))

    if missing_list:

        pd.Series(
            sorted(set(missing_list))
        ).to_csv(
            EXTRACT_DIR /
            f"{bed_file.stem}_missing_ids.txt",
            index=False,
            header=False
        )

    print(f"Extracted LCRs : {len(out_df)}")
    print(f"Missing IDs    : {missing_ids}")
    print(f"Invalid coords : {invalid_coords}")
    print("start<1 :", bad_start)
    print("end>len :", bad_end)
    print("start>end :", bad_order)

print("[SUCCESS] Run_pipeline_1 completed")
print("\nPipeline step 1 complete.")