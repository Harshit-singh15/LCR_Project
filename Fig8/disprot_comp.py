import pandas as pd
from collections import Counter

# =====================================================
# INPUT / OUTPUT
# =====================================================

input_file = r"Fig8\DisProt release_2025_12 with_ambiguous_evidences.tsv"

output_file = r"fig8A_disprot_complexity.tsv"

# =====================================================


def find_most_frequent_aa(sequence):

    if len(sequence) == 0:
        return "", 0

    freqs = Counter(sequence)

    aa, count = freqs.most_common(1)[0]

    return aa, (count / len(sequence)) * 100


def find_kmers(sequence):

    kmer_counts = {}

    for k in range(1, len(sequence)):

        seen = Counter()

        for i in range(len(sequence) - k + 1):

            kmer = sequence[i:i+k]

            seen[kmer] += 1

        kmer_counts[k] = {
            kmer: count
            for kmer, count in seen.items()
            if count > 1
        }

    return kmer_counts


def min_mutation_percent(sequence, kmer_counts):

    min_mutations = len(sequence)

    best_kmer = None

    for k, kmers in kmer_counts.items():

        for kmer in kmers:

            repetitions = sequence.count(kmer)

            required_mutations = (
                len(sequence)
                - (repetitions * len(kmer))
            )

            if required_mutations < min_mutations:

                min_mutations = required_mutations

                best_kmer = kmer

    if best_kmer:

        return (
            min_mutations / len(sequence)
        ) * 100

    return 100


print("Loading DisProt file...")

df = pd.read_csv(
    input_file,
    sep="\t"
)

print(f"Regions: {len(df):,}")

results = []

for idx, row in df.iterrows():

    sequence = str(row["region_sequence"])

    if sequence == "nan":
        continue

    _, most_freq_percent = find_most_frequent_aa(
        sequence
    )

    kmers = find_kmers(sequence)

    mutation_percent = min_mutation_percent(
        sequence,
        kmers
    )

    results.append([
        row["acc"],
        row["region_id"],
        row["ec_name"],
        mutation_percent,
        most_freq_percent
    ])

    if idx % 100 == 0:
        print(idx)

out = pd.DataFrame(
    results,
    columns=[
        "acc",
        "region_id",
        "ec_name",
        "Mutation_Percent",
        "Most_Frequent_AA_Percent"
    ]
)

out.to_csv(
    output_file,
    sep="\t",
    index=False
)

print("\nSaved:")
print(output_file)