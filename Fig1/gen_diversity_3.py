from pathlib import Path
from collections import Counter
import pandas as pd
import math
import gc
import sys

INPUT_DIR = sys.argv[1]      # Extracted Sequences
OUTPUT_DIR = sys.argv[2]     # Shannon entropy output

Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)


def shannon_entropy(seq):

    if not seq:
        return 0.0

    counts = Counter(seq)
    seq_len = len(seq)

    entropy = 0.0

    for count in counts.values():

        p = count / seq_len

        entropy -= p * math.log2(p)

    return entropy


for file in Path(INPUT_DIR).glob("*_lcrs.tsv"):

    print(f"\nProcessing {file.name}")

    df = pd.read_csv(
        file,
        sep="\t",
        usecols=["Sequence"],
        dtype={"Sequence": "string"}
    )

    out_file = (
        Path(OUTPUT_DIR)
        /
        f"{file.stem.replace('_lcrs','')}_SNS.tsv"
    )

    total = 0.0
    maximum = float("-inf")
    count = 0

    with open(out_file, "w") as out:

        out.write("Shannon_Entropy\n")

        for seq in df["Sequence"]:

            entropy = shannon_entropy(str(seq))

            out.write(f"{entropy}\n")

            total += entropy

            if entropy > maximum:
                maximum = entropy

            count += 1

    mean_entropy = total / count if count else 0

    print(f"Mean = {mean_entropy:.3f}")
    print(f"Max = {maximum:.3f}")

    del df
    gc.collect()

print("\nDone.")