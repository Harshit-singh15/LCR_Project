# infer_missing_residues.py

from collections import defaultdict

# -------------------------
# protein lengths
# -------------------------

lengths = {}

with open("Fig8\mouse_protein_lengths.sorted.tsv") as f:
    for line in f:
        acc, length = line.rstrip().split("\t")
        lengths[acc] = int(length)

# -------------------------
# observed intervals
# -------------------------

observed = defaultdict(list)

with open("Fig8\mouse_sifts_observed.sorted.tsv") as f:
    for line in f:

        cols = line.rstrip().split("\t")

        acc = cols[2]

        sp_beg = int(cols[7])
        sp_end = int(cols[8])

        observed[acc].append((sp_beg, sp_end))

# -------------------------
# missing intervals
# -------------------------

out = open("Fig8\\mouse_missing_residues.bed", "w")

for acc in observed:

    intervals = sorted(observed[acc])

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

    for i in range(len(merged)-1):

        gap_start = merged[i][1] + 1
        gap_end   = merged[i+1][0] - 1

        if gap_end >= gap_start:

            gap_len = gap_end - gap_start + 1

            if gap_len >= 10:

                out.write(
                    f"{acc}\t{gap_start}\t{gap_end}\n"
                )
out.close()

print("done")