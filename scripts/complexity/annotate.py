import pandas as pd

# ===== EDIT THESE =====

input_file = r"reference_lcr\mouse_windows_classified.tsv"

# ======================


def resolve_overlaps():

    print("Reading classified windows...")

    df = pd.read_csv(
        input_file,
        sep="\t"
    )

    position_annotations = {}

    print("Expanding windows to positions...")

    for _, row in df.iterrows():

        protein = row["Protein_ID"]

        start = int(row["Start"])

        end = int(row["End"])

        classification = row["Classification"]

        for pos in range(start, end + 1):

            key = (protein, pos)

            if key not in position_annotations:
                position_annotations[key] = set()

            position_annotations[key].add(
                classification
            )

    print("Resolving overlaps...")

    resolved_annotations = {}

    for key, classifications in position_annotations.items():

        if "LCR" in classifications:

            resolved_annotations[key] = "LCR"

        elif "CBR" in classifications:

            resolved_annotations[key] = "CBR"

        else:

            resolved_annotations[key] = "HCR"

    rows = []

    for (protein, pos), classification in resolved_annotations.items():

        rows.append(
            [
                protein,
                pos,
                classification
            ]
        )

    resolved_df = pd.DataFrame(
        rows,
        columns=[
            "Protein_ID",
            "Position",
            "Classification"
        ]
    )

    resolved_df = resolved_df.sort_values(
        ["Protein_ID", "Position"]
    )

    resolved_df.to_csv(
        input_file,
        sep="\t",
        index=False
    )

    print(f"\nUpdated: {input_file}")

    print(
        f"Total positions: {len(resolved_df):,}"
    )


if __name__ == "__main__":
    resolve_overlaps()