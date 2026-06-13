from Bio import SeqIO
from pathlib import Path
import pandas as pd

# =====================================================
# CONFIG
# =====================================================

FASTA_FILE = "outputs\\mouse_proteome.fasta"

BED_FOLDER = "bed_files"

OUTPUT_DIR = "outputs_prerequisite"

# =====================================================
# CREATE OUTPUT DIRECTORIES
# =====================================================

Path(OUTPUT_DIR).mkdir(exist_ok=True)

EXTRACT_DIR = Path(OUTPUT_DIR) / "extracted_sequences"
EXTRACT_DIR.mkdir(exist_ok=True)

# =====================================================
# LOAD FASTA
# =====================================================

print("\nLoading proteome...")

proteins = {}
protein_lengths = {}

for record in SeqIO.parse(FASTA_FILE, "fasta"):

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

            mouse_name = parts[2]

            proteins[accession] = sequence
            proteins[mouse_name] = sequence

            protein_lengths[accession] = length
            protein_lengths[mouse_name] = length

print(f"Loaded {len(protein_lengths)} identifiers")

# =====================================================
# SAVE PROTEIN LENGTHS
# =====================================================

canonical_lengths = []

for record in SeqIO.parse(FASTA_FILE, "fasta"):

    full_id = record.description.split()[0]

    canonical_lengths.append(
        [full_id, len(record.seq)]
    )

length_df = pd.DataFrame(
    canonical_lengths,
    columns=[
        "Protein_ID",
        "Protein_Length"
    ]
)

length_df.to_csv(
    Path(OUTPUT_DIR) / "protein_lengths.tsv",
    sep="\t",
    index=False
)

print("Saved protein_lengths.tsv")

# =====================================================
# PROCESS BED FILES
# =====================================================

bed_files = sorted(Path(BED_FOLDER).glob("*.bed"))

print(f"\nFound {len(bed_files)} BED files")

for bed_file in bed_files:

    print(f"\nProcessing {bed_file.name}")

    try:

        df = pd.read_csv(
            bed_file,
            sep="\t"
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
            continue

        protein_seq = proteins[protein_id]

        protein_len = len(protein_seq)

        # ----------------------------------------
        # coordinate validation
        # ----------------------------------------

        if start < 1:

            invalid_coords += 1
            continue

        if end > protein_len:

            invalid_coords += 1
            continue

        if start > end:

            invalid_coords += 1
            continue

        # ----------------------------------------
        # 1-based inclusive extraction
        # ----------------------------------------

        lcr_seq = protein_seq[start - 1:end]

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

    output_file = (
        EXTRACT_DIR /
        f"{bed_file.stem}_lcrs.tsv"
    )

    out_df.to_csv(
        output_file,
        sep="\t",
        index=False
    )

    print(f"Extracted LCRs : {len(out_df)}")
    print(f"Missing IDs    : {missing_ids}")
    print(f"Invalid coords : {invalid_coords}")

print("\nPipeline step 1 complete.")