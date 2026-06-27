import math
from collections import Counter
from Bio import SeqIO

# ===== EDIT THESE =====

fasta_file = r"celegans\celegans_cleaned.fasta"

bed_file = r"celegans\dataforFig6\celegans_windows.bed"

output_file = r"celegans\dataforFig6\celegans_windows_out.tsv"

# ======================


def read_bed(file):

    bed_regions = []

    with open(file, "r") as f:

        for line in f:

            parts = line.strip().split()

            if len(parts) >= 3:

                bed_regions.append(
                    (
                        parts[0],
                        int(parts[1]),
                        int(parts[2])
                    )
                )

    return bed_regions


def extract_sequences(fasta_file, bed_regions):

    sequences = {}

    for record in SeqIO.parse(fasta_file, "fasta"):

        sequences[record.id] = str(record.seq)

    extracted = []

    for protein_id, start, end in bed_regions:

        if protein_id not in sequences:
            continue

        seq = sequences[protein_id]

        # BED is 1-based inclusive
        subseq = seq[start - 1:end]

        if len(subseq) == 0:
            continue

        extracted.append(
            (
                protein_id,
                start,
                end,
                subseq
            )
        )

    return extracted


def calculate_entropy(sequence):

    if len(sequence) == 0:
        return 0.0

    length = len(sequence)

    freqs = Counter(sequence)

    return -sum(
        (count / length) * math.log2(count / length)
        for count in freqs.values()
    )


def find_most_frequent_aa(sequence):

    if len(sequence) == 0:
        return "NA", 0.0

    freqs = Counter(sequence)

    most_common_aa, most_common_count = freqs.most_common(1)[0]

    return (
        most_common_aa,
        (most_common_count / len(sequence)) * 100
    )


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
        return None, 100.0, ""

    min_mutations = len(sequence)

    best_kmer = None

    mutation_data = []

    for k, kmers in kmer_counts.items():

        for kmer in kmers:

            repetitions = sequence.count(kmer)

            required_mutations = (
                len(sequence)
                - (repetitions * len(kmer))
            )

            mutation_percent = (
                required_mutations / len(sequence)
            ) * 100

            mutation_data.append(
                f"{kmer}:{mutation_percent:.2f}"
            )

            if required_mutations < min_mutations:

                min_mutations = required_mutations

                best_kmer = kmer

    overall_mutation_percent = (
        (min_mutations / len(sequence)) * 100
        if best_kmer
        else 100.0
    )

    return (
        best_kmer,
        overall_mutation_percent,
        ",".join(mutation_data)
    )


def main():

    print("Reading BED file...")
    bed_regions = read_bed(bed_file)

    print("Extracting sequences...")
    extracted_sequences = extract_sequences(
        fasta_file,
        bed_regions
    )

    print(f"Processing {len(extracted_sequences)} windows...")

    with open(output_file, "w") as out:

        out.write(
            "Protein_ID\tStart\tEnd\tEntropy\tMost_Frequent_AA\tMost_Frequent_AA_Percent\tKmers_At_Least_Twice\tBest_Kmer\tMutation_Percent\tKmer_Mutation_List\tSequence\n"
        )

        for (
            protein_id,
            start,
            end,
            sequence
        ) in extracted_sequences:

            entropy = calculate_entropy(sequence)

            (
                most_freq_aa,
                most_freq_aa_percent
            ) = find_most_frequent_aa(sequence)

            kmer_counts = find_kmers(sequence)

            (
                best_kmer,
                mutation_percent,
                mutation_list
            ) = min_mutation_percent(
                sequence,
                kmer_counts
            )

            kmer_list = ",".join(
                [
                    kmer
                    for kmers in kmer_counts.values()
                    for kmer in kmers
                ]
            )

            out.write(
                f"{protein_id}\t"
                f"{start}\t"
                f"{end}\t"
                f"{entropy:.4f}\t"
                f"{most_freq_aa}\t"
                f"{most_freq_aa_percent:.2f}\t"
                f"{kmer_list}\t"
                f"{best_kmer}\t"
                f"{mutation_percent:.2f}\t"
                f"{mutation_list}\t"
                f"{sequence}\n"
            )

    print(f"\nSaved: {output_file}")


if __name__ == "__main__":
    main()