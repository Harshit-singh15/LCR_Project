from pathlib import Path
import pandas as pd
import re

# ==========================================================
# Directories
# ==========================================================

INPUT_DIR = Path(
    r"fruitfly\dataforFig2\06_peptide_counts"
)

OUTPUT_DIR = Path(
    r"fruitfly\dataforFig2\07_motif_list"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

# ==========================================================
# Read all peptide count files
# ==========================================================

motif_matrix = {}

consensus_levels = []

for file in sorted(
    INPUT_DIR.glob("*_peptide_counts.tsv")
):

    match = re.search(
        r"consensus_(\d+)",
        file.stem
    )

    if match is None:
        continue

    consensus = int(match.group(1))

    consensus_levels.append(consensus)

    print(f"Processing Consensus {consensus}")

    df = pd.read_csv(
        file,
        sep="\t"
    )

    for _, row in df.iterrows():

        motif = row["Best-Peptide"]

        count = row["Count"]

        if motif not in motif_matrix:

            motif_matrix[motif] = {}

        motif_matrix[motif][
            consensus
        ] = count

# ==========================================================
# Convert to DataFrame
# ==========================================================

consensus_levels = sorted(
    consensus_levels
)

matrix = pd.DataFrame.from_dict(
    motif_matrix,
    orient="index"
)

matrix = matrix.reindex(
    columns=consensus_levels
)

matrix = matrix.fillna(0)

matrix = matrix.astype(int)

matrix.index.name = "Motif"

matrix.columns = [
    f"Consensus_{i}"
    for i in consensus_levels
]

# ==========================================================
# Sort motifs
# ==========================================================

matrix["Total"] = matrix.sum(axis=1)

matrix = matrix.sort_values(
    "Total",
    ascending=False
)

# ==========================================================
# Save
# ==========================================================

outfile = OUTPUT_DIR / "motifs.tsv"

matrix.to_csv(
    outfile,
    sep="\t"
)

print("\n================================")

print(
    f"Unique motifs : {len(matrix)}"
)

print(
    f"Saved : {outfile}"
)

print("================================")