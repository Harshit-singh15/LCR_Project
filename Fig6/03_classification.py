from pathlib import Path
import pandas as pd
import sys

input_file = Path(sys.argv[1])  # Input file with metrics
output_file = Path(sys.argv[2])  # Output file for classified regions

output_file.parent.mkdir(
    parents=True,
    exist_ok=True
)

def classify_region(most_freq_aa_percent, mutation_percent):
    if most_freq_aa_percent < 50 and mutation_percent <= 50:
        return "CBR"
    elif most_freq_aa_percent >= 50 and mutation_percent <= 50:
        return "LCR"
    elif most_freq_aa_percent < 50 and mutation_percent > 50:
        return "HCR"
    else:
        return "Unknown"

def main():
    print("Reading and processing file in chunks...")

    # Define columns to keep
    columns_to_keep = [
        "Protein_ID",
        "Start",
        "End",
        "Most_Frequent_AA_Percent",
        "Mutation_Percent"
    ]
    
    # Track overall classification counts across all chunks
    from collections import Counter
    total_counts = Counter()

    # Flag to handle writing headers on the first chunk only
    is_first_chunk = True

    # Adjust chunksize (e.g., 50,000 or 100,000 rows) depending on your needs
    chunk_size = 50000 

    # Read the file iteratively
    chunks = pd.read_csv(
        input_file,
        sep="\t",
        usecols=columns_to_keep,  # Memory saver: only load required columns
        chunksize=chunk_size
    )

    for chunk in chunks:
        # Classify the current chunk
        chunk["Classification"] = chunk.apply(
            lambda row: classify_region(
                row["Most_Frequent_AA_Percent"],
                row["Mutation_Percent"]
            ),
            axis=1
        )

        # Update the aggregate counts
        total_counts.update(chunk["Classification"])

        # Append to the output file
        chunk.to_csv(
            output_file,
            sep="\t",
            index=False,
            mode='a',                  # 'a' stands for append mode
            header=is_first_chunk      # Only write the header for the very first chunk
        )
        
        is_first_chunk = False

    print(f"\nSaved: {output_file}")
    print("\nClassification counts:")
    for classification, count in total_counts.items():
        print(f"{classification:<10} {count}")

if __name__ == "__main__":
    main()