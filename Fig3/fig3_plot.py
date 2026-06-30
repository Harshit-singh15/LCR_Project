import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

INPUT_DIR = Path(
    "Fig3/02_metrics"
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

    method = (
        file.stem
        .replace("_purity","")
    )

    df = pd.read_csv(
        file,
        sep="\t"
    )

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
    "Purity distribution across LCR detection methods"
)

plt.legend(
    bbox_to_anchor=(1.02,1),
    loc="upper left"
)

plt.grid(
    alpha=0.3
)

plt.tight_layout()

plt.savefig(
    "Fig3/03_plots/Fig3_purity_distribution.png",
    dpi=600,
    bbox_inches="tight"
)

plt.show()