import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import sys

INPUT_DIR = Path(
    sys.argv[1]
)

thresholds = np.arange(
    0,
    1.01,
    0.01
)

plt.figure(
    figsize=(10,7)
)

for file in INPUT_DIR.glob(
    "*_purity.tsv"
):

    method = file.stem.replace("_purity", "")

# Remove organism suffix
    for organism in [
        "_celegans",
        "_mouse",
        "_zebrafish",
        "_human",
        "_yeast",
        "_Fruitfly",
        "_arabidopsis",
        "_ecoli"
    ]:
        method = method.replace(organism, "")

    if file.stat().st_size == 0:
        print(f"Skipping empty file: {file.name}")
        continue

    try:
        df = pd.read_csv(file, sep="\t")
    except pd.errors.EmptyDataError:
        print(f"Skipping empty file: {file.name}")
        continue

    purity = (
        df["Purity"]
        .values
    )

    retained = []

    total = len(purity)

    for t in thresholds:

        retained.append(
            np.sum(
                purity >= t
            ) / total
        )

    plt.plot(
        thresholds,
        retained,
        linewidth=2,
        label=method
    )

plt.xlabel(
    "Purity threshold"
)

plt.ylabel(
    "Proportion of LCRs retained"
)

plt.title(
    "Fig 3: E. coli : Purity distribution across LCR detection methods"
)

plt.legend(
    bbox_to_anchor=(1.02,1),
    loc="upper left"
)

plt.grid(
    alpha=0.3
)

plt.tight_layout()

output_path = Path(sys.argv[2])
output_path.parent.mkdir(parents=True, exist_ok=True)

plt.savefig(
    output_path/"Fig3_Purity.png",
    dpi=600,
    bbox_inches="tight"
)
