import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pathlib import Path
import sys

INPUT_DIR = Path(sys.argv[1])
OUTPUT_DIR = Path(sys.argv[2])
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

organisms = [
    "_mouse","_zebrafish","_human","_yeast",
    "_celegans","_Fruitfly","_arabidopsis","_ecoli"
]

def clean_name(name):
    name = name.replace("_metrics", "")
    for org in organisms:
        name = name.replace(org, "")
    return name

for file in sorted(INPUT_DIR.glob("*_metrics.tsv")):

    method = clean_name(file.stem)
    print(f"Processing: {file.name}")

    try:
        df = pd.read_csv(file, sep="\t")

    except Exception as e:
        print(f"[SKIPPED] {file.name}")
        print(f"Reason : {type(e).__name__}: {e}")

        try:
            print(f"Size : {file.stat().st_size} bytes")
            with open(file, "rb") as f:
                print(f"First bytes : {repr(f.read(80))}")
        except:
            pass

        print()
        continue

    df["MutBin"] = (df["Mutation_Percent"] / 10).round().astype(int).clip(0, 10)
    df["FreqBin"] = (df["Most_Common_AA_Percent"] / 10).round().astype(int).clip(0, 10)

    heatmap = np.zeros((11, 11))

    for _, row in df.iterrows():
        heatmap[int(row["FreqBin"]), int(row["MutBin"])] += 1

    fig, ax = plt.subplots(figsize=(7, 6))

    im = ax.imshow(
        heatmap,
        origin="lower",
        cmap="Reds",
        norm=LogNorm(vmin=1, vmax=max(1, heatmap.max())),
        aspect="auto"
    )

    ax.set_xticks(range(11))
    ax.set_yticks(range(11))
    ax.set_xticklabels([i * 10 for i in range(11)], rotation=45)
    ax.set_yticklabels([i * 10 for i in range(11)])

    ax.set_xlabel("Mutation Percentage")
    ax.set_ylabel("Most Frequent AA Percentage")
    ax.set_title(method)

    plt.colorbar(im, ax=ax, label="Number of LCRs")
    plt.tight_layout()

    plt.savefig(
        OUTPUT_DIR / f"{method}_Fig5.png",
        dpi=600,
        bbox_inches="tight"
    )

    plt.close()

print("Done.")