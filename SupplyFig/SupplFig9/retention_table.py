import os
import pandas as pd
import numpy as np

# =====================================
# EDIT THESE
# =====================================

INPUT_DIR = r"SupplyFig\SupplFig9\complexity_xstream"

OUTPUT_FILE = r"SupplyFig\SupplFig9\retention_table.tsv"

# =====================================

thresholds = np.arange(0, 100, 10)

results = []

tsv_files = sorted([
    f for f in os.listdir(INPUT_DIR)
    if f.endswith(".tsv")
])

print("Files found:", len(tsv_files))

for file in tsv_files:

    path = os.path.join(
        INPUT_DIR,
        file
    )

    tool = os.path.splitext(file)[0]

    print("Processing:", tool)

    df = pd.read_csv(
        path,
        sep="\t"
    )

    if "Most_Frequent_AA_Percent" not in df.columns:

        print(
            "Skipping:",
            file,
            "(column missing)"
        )

        continue

    total = len(df)

    if total == 0:

        print(
            "Skipping:",
            file,
            "(empty file)"
        )

        continue

    for purity in thresholds:

        retained = (
            df["Most_Frequent_AA_Percent"]
            >= purity
        ).sum()

        results.append([
            tool,
            purity,
            retained,
            total,
            retained / total
        ])

out = pd.DataFrame(
    results,
    columns=[
        "Tool",
        "Purity",
        "Retained",
        "Total",
        "Proportion"
    ]
)

out.to_csv(
    OUTPUT_FILE,
    sep="\t",
    index=False
)

print("\nSaved:", OUTPUT_FILE)
print(out.head())