import pandas as pd
from collections import Counter

# =====================================================
# INPUT / OUTPUT
# =====================================================

input_file = r"zebrafish\dataforFig8\Disprot_zebrafish.tsv"

output_file = r"zebrafish\dataforFig8\fig8A_disprot_complexity.tsv"

log_file = r"zebrafish\dataforFig8\fig8A_disprot_complexity.log"

# =====================================================


def clean_sequence(sequence):
    """
    Remove ambiguous amino acids.

    This is a defensive preprocessing step.
    Standard UniProt sequences are expected to be
    unaffected.
    """

    allowed = set("ACDEFGHIKLMNPQRSTVWY")

    return "".join(
        aa
        for aa in sequence.upper()
        if aa in allowed
    )


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

print(f"Regions : {len(df):,}")

results = []

processed = 0
skipped = 0
cleaned = 0
errors = 0

log = open(log_file, "w")

for index, row in df.iterrows():

    try:

        sequence = str(row["Region sequence"])

        if sequence == "nan":

            skipped += 1

            log.write(
                f"Skipped missing sequence : "
                f"{row['Region ID']}\n"
            )

            continue

        original_length = len(sequence)

        sequence = clean_sequence(sequence)

        if len(sequence) != original_length:
            cleaned += 1

        _, most_freq_percent = (
            find_most_frequent_aa(sequence)
        )

        kmers = find_kmers(sequence)

        mutation_percent = (
            min_mutation_percent(
                sequence,
                kmers
            )
        )

        results.append([

            row["UniProt ACC"],

            row["Region ID"],

            row["ECO Term name"],

            mutation_percent,

            most_freq_percent

        ])

        processed += 1

        if processed % 100 == 0:
            print(processed)

    except Exception as e:

        errors += 1

        log.write(
            f"ERROR : "
            f"{row['Region ID']} : "
            f"{e}\n"
        )

log.close()

out = pd.DataFrame(

    results,

    columns=[

        "Protein_ID",

        "Region_ID",

        "Experimental_Method",

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

print("\nStatistics")

print("------------------------")

print(f"Total regions     : {len(df):,}")

print(f"Processed         : {processed:,}")

print(f"Skipped           : {skipped:,}")

print(f"Sequences cleaned : {cleaned:,}")

print(f"Errors            : {errors:,}")