import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pathlib import Path
import sys

# =====================================================
# INPUT / OUTPUT
# =====================================================

INPUT_DIR = Path(sys.argv[1])
OUTPUT_DIR = Path(sys.argv[2])

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

# =====================================================
# Helper
# =====================================================

def clean_name(name):
    """
    Remove filename suffix used by the pipeline.
    """
    return name.replace("_metrics", "")

# =====================================================
# Process every metrics file
# =====================================================

files = sorted(INPUT_DIR.glob("*_metrics.tsv"))

if not files:
    raise FileNotFoundError(
        f"No *_metrics.tsv files found in {INPUT_DIR}"
    )

for file in files:

    method = clean_name(file.stem)

    print(f"Processing: {file.name}")

    try:
        df = pd.read_csv(
            file,
            sep="\t"
        )

    except Exception as e:

        print(f"[SKIPPED] {file.name}")
        print(f"Reason: {type(e).__name__}: {e}")

        try:
            print(f"Size: {file.stat().st_size} bytes")
            with open(file, "rb") as f:
                print(f"First bytes: {repr(f.read(80))}")
        except Exception:
            pass

        print()
        continue

    required = {
        "Mutation_Percent",
        "Most_Common_AA_Percent"
    }

    if not required.issubset(df.columns):
        print(
            f"Skipping {file.name}: Missing columns "
            f"{required - set(df.columns)}"
        )
        continue

    df["MutBin"] = (
        df["Mutation_Percent"] / 10
    ).round().astype(int).clip(0, 10)

    df["FreqBin"] = (
        df["Most_Common_AA_Percent"] / 10
    ).round().astype(int).clip(0, 10)

    heatmap = np.zeros((11, 11))

    for _, row in df.iterrows():
        heatmap[
            int(row["FreqBin"]),
            int(row["MutBin"])
        ] += 1

    fig, ax = plt.subplots(
        figsize=(7, 6)
    )

    im = ax.imshow(
        heatmap,
        origin="lower",
        cmap="Reds",
        norm=LogNorm(
            vmin=1,
            vmax=max(1, heatmap.max())
        ),
        aspect="auto"
    )

    ax.set_xticks(range(11))
    ax.set_yticks(range(11))

    ax.set_xticklabels(
        [i * 10 for i in range(11)],
        rotation=45
    )

    ax.set_yticklabels(
        [i * 10 for i in range(11)]
    )

    ax.set_xlabel(
        "Mutation Percentage"
    )

    ax.set_ylabel(
        "Most Frequent AA Percentage"
    )

    # Uses filename as panel title
    ax.set_title(method)

    plt.colorbar(
        im,
        ax=ax,
        label="Number of LCRs"
    )

    plt.tight_layout()

    output_file = OUTPUT_DIR / f"{method}_Fig5.png"

    plt.savefig(
        output_file,
        dpi=600,
        bbox_inches="tight"
    )

    plt.close(fig)

print("Done.")