from Bio import SeqIO
from collections import Counter
import pandas as pd

# =====================================================
# INPUT / OUTPUT
# =====================================================

input_fasta = r"celegans\dataforFig8\celegans_missing_residues.fa"

output_file = r"celegans\dataforFig8\fig8B_missing_residue_complexity.tsv"

log_file = r"celegans\dataforFig8\fig8B_missing_complexity.log"

# =====================================================


def clean_sequence(sequence):
    """
    Remove ambiguous amino acids.

    Defensive preprocessing.
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

    if len(sequence) < 2:
        return {}

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

    if len(sequence) == 0:
        return 100

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


print("Reading FASTA...")

results = []

processed = 0
cleaned = 0
errors = 0

log = open(log_file, "w")

for record in SeqIO.parse(input_fasta, "fasta"):

    try:

        sequence = str(record.seq)

        original_length = len(sequence)

        sequence = clean_sequence(sequence)

        if len(sequence) != original_length:

            cleaned += 1

        acc, start, end = record.id.split("|")

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

            acc,

            int(start),

            int(end),

            len(sequence),

            mutation_percent,

            most_freq_percent

        ])

        processed += 1

        if processed % 100 == 0:

            print(processed)

    except Exception as e:

        errors += 1

        log.write(
            f"{record.id}\t{e}\n"
        )

log.close()

out = pd.DataFrame(

    results,

    columns=[

        "Protein_ID",

        "Start",

        "End",

        "Length",

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

print("-----------------------")

print(f"Regions processed : {processed:,}")

print(f"Sequences cleaned : {cleaned:,}")

print(f"Errors            : {errors:,}")