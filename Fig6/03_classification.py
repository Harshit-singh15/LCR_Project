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

    print("Reading file...")

    df = pd.read_csv(
        input_file,
        sep="\t"
    )

    print("Classifying windows...")

    df["Classification"] = df.apply(
        lambda row: classify_region(
            row["Most_Frequent_AA_Percent"],
            row["Mutation_Percent"]
        ),
        axis=1
    )

    output_df = df[
        [
            "Protein_ID",
            "Start",
            "End",
            "Most_Frequent_AA_Percent",
            "Mutation_Percent",
            "Classification"
        ]
    ]

    output_df.to_csv(
        output_file,
        sep="\t",
        index=False
    )

    print(f"\nSaved: {output_file}")

    print("\nClassification counts:")

    print(
        output_df["Classification"]
        .value_counts()
    )


if __name__ == "__main__":
    main()