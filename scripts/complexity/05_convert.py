import pandas as pd

# ===== EDIT THESE =====

input_file = r"celegans\dataforFig6\celegans_windows_classified.tsv"

output_file = r"celegans\dataforFig6\celegans_windows_real.bed"

# ======================


def find_contiguous_regions(df):

    results = []

    for protein, group in df.groupby(
        "Protein_ID",
        sort=False,
        observed=False
    ):

        start = None
        current_class = None

        for _, row in group.iterrows():

            pos = row["Position"]
            classification = row["Classification"]

            if start is None:

                start = pos
                current_class = classification

            elif classification != current_class:

                results.append(
                    [
                        protein,
                        start,
                        prev_pos,
                        current_class
                    ]
                )

                start = pos
                current_class = classification

            prev_pos = pos

        results.append(
            [
                protein,
                start,
                prev_pos,
                current_class
            ]
        )

    return results


def main():

    print("Reading position annotations...")

    df = pd.read_csv(
        input_file,
        sep="\t"
    )

    df = df[
        df["Position"] > 0
    ]

    df["Protein_ID"] = pd.Categorical(
        df["Protein_ID"],
        categories=df["Protein_ID"].unique(),
        ordered=True
    )

    print("Building contiguous regions...")

    regions = find_contiguous_regions(df)

    output_df = pd.DataFrame(
        regions,
        columns=[
            "Protein_ID",
            "Start_Position",
            "End_Position",
            "Classification"
        ]
    )

    output_df.to_csv(
        output_file,
        sep="\t",
        index=False
    )

    print(f"\nSaved: {output_file}")

    print(
        f"Total regions: {len(output_df):,}"
    )

    print("\nClassification counts:")

    print(
        output_df["Classification"]
        .value_counts()
    )


if __name__ == "__main__":
    main()