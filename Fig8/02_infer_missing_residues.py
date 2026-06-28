from collections import defaultdict
import pandas as pd

# ==========================================================
# INPUT FILES
# ==========================================================

protein_lengths = r"celegans\dataforFig1\protein_lengths.tsv"

observed_file = r"celegans\dataforFig8\celegans_sifts_observed.tsv"

output_bed = r"celegans\dataforFig8\celegans_missing_residues.bed"

minimum_gap = 10

# ==========================================================


print("Loading protein lengths...")

lengths = {}

with open(protein_lengths) as f:
    next(f)  # skip header
    for line in f:

        protein, length = line.rstrip().split("\t")

        lengths[protein] = int(length)

print(f"Proteins loaded : {len(lengths):,}")

print("\nLoading observed residue intervals...")

df = pd.read_csv(
    observed_file,
    sep="\t"
)

observed = defaultdict(list)

missing_lengths = 0

for _, row in df.iterrows():

    acc = row["SP_PRIMARY"]

    start = int(row["SP_BEG"])

    end = int(row["SP_END"])

    if acc not in lengths:
        continue

    if start > end:
        continue

    observed[acc].append((start, end))

print(f"Proteins with observations : {len(observed):,}")

print("\nInferring missing regions...")

out = open(output_bed, "w")

gap_count = 0

for acc, intervals in observed.items():

    protein_length = lengths[acc]

    intervals = sorted(intervals)

    merged = []

    for start, end in intervals:

        if not merged:

            merged.append([start, end])

        elif start <= merged[-1][1] + 1:

            merged[-1][1] = max(
                merged[-1][1],
                end
            )

        else:

            merged.append([start, end])

    # -----------------------------------------
    # N-terminal gap
    # -----------------------------------------

    if merged[0][0] > 1:

        gap_start = 1

        gap_end = merged[0][0] - 1

        if gap_end - gap_start + 1 >= minimum_gap:

            out.write(
                f"{acc}\t{gap_start}\t{gap_end}\n"
            )

            gap_count += 1

    # -----------------------------------------
    # Internal gaps
    # -----------------------------------------

    for i in range(len(merged) - 1):

        gap_start = merged[i][1] + 1

        gap_end = merged[i + 1][0] - 1

        if gap_end >= gap_start:

            if gap_end - gap_start + 1 >= minimum_gap:

                out.write(
                    f"{acc}\t{gap_start}\t{gap_end}\n"
                )

                gap_count += 1

    # -----------------------------------------
    # C-terminal gap
    # -----------------------------------------

    if merged[-1][1] < protein_length:

        gap_start = merged[-1][1] + 1

        gap_end = protein_length

        if gap_end - gap_start + 1 >= minimum_gap:

            out.write(
                f"{acc}\t{gap_start}\t{gap_end}\n"
            )

            gap_count += 1

out.close()

print("\nDone.")

print(f"Missing regions : {gap_count:,}")

print("\nSaved:")

print(output_bed)