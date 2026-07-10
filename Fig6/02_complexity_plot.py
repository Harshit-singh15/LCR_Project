import math
from collections import Counter
from pathlib import Path
from Bio import SeqIO
import sys

# ===== EDIT THESE =====
fasta_file = Path(sys.argv[1])  # Proteome FASTA file
bed_file = Path(sys.argv[2])   # Input BED file
output_file = Path(sys.argv[3]) # Output file for metrics
# ======================

output_file.parent.mkdir(
    parents=True,
    exist_ok=True
)

def read_bed_generator(file):
    """Yields BED regions one line at a time to save memory."""
    with open(file, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                yield parts[0], int(parts[1]), int(parts[2])

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
    return most_common_aa, (most_common_count / len(sequence)) * 100

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
            required_mutations = len(sequence) - (repetitions * len(kmer))
            mutation_percent = (required_mutations / len(sequence)) * 100
            mutation_data.append(f"{kmer}:{mutation_percent:.2f}")

            if required_mutations < min_mutations:
                min_mutations = required_mutations
                best_kmer = kmer

    overall_mutation_percent = (
        (min_mutations / len(sequence)) * 100
        if best_kmer
        else 100.0
    )
    return best_kmer, overall_mutation_percent, ",".join(mutation_data)

def main():
    print("Indexing FASTA file (Memory efficient)...")
    # SeqIO.index acts like a dictionary but keeps sequences on disk, loading them on-demand
    fasta_index = SeqIO.index(str(fasta_file), "fasta")

    print("Processing windows and writing output on the fly...")
    
    with open(output_file, "w") as out:
        out.write(
            "Protein_ID\tStart\tEnd\tEntropy\tMost_Frequent_AA\tMost_Frequent_AA_Percent\t"
            "Kmers_At_Least_Twice\tBest_Kmer\tMutation_Percent\tKmer_Mutation_List\tSequence\n"
        )

        processed_count = 0

        # Stream regions line by line from BED
        for protein_id, start, end in read_bed_generator(bed_file):
            if protein_id not in fasta_index:
                continue

            # Load ONLY the record we need from disk into memory
            record = fasta_index[protein_id]
            seq = str(record.seq)

            # BED is 1-based inclusive
            subseq = seq[start - 1:end]
            if len(subseq) == 0:
                continue

            # Run metrics
            entropy = calculate_entropy(subseq)
            most_freq_aa, most_freq_aa_percent = find_most_frequent_aa(subseq)
            kmer_counts = find_kmers(subseq)
            best_kmer, mutation_percent, mutation_list = min_mutation_percent(subseq, kmer_counts)

            kmer_list = ",".join(
                [kmer for kmers in kmer_counts.values() for kmer in kmers]
            )

            # Write immediately to file
            out.write(
                f"{protein_id}\t{start}\t{end}\t{entropy:.4f}\t"
                f"{most_freq_aa}\t{most_freq_aa_percent:.2f}\t"
                f"{kmer_list}\t{best_kmer}\t{mutation_percent:.2f}\t"
                f"{mutation_list}\t{subseq}\n"
            )

            processed_count += 1
            if processed_count % 1000 == 0:
                print(f"Processed {processed_count} regions...")

    print(f"\nSaved: {output_file}")

if __name__ == "__main__":
    main()