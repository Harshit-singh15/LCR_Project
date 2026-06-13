from pathlib import Path
import pandas as pd

TOTAL_PROTEINS = 21853

INPUT_DIR = Path("outputs_prerequisite\\LCR_Count")

print("\n=== COVERAGE QC ===\n")

for file in INPUT_DIR.glob("*.tsv"):

    df = pd.read_csv(file, sep="\t")

    total = df["Count"].sum()

    status = "PASS"

    if total != TOTAL_PROTEINS:
        status = "FAIL"

    print(
        f"{file.stem:40s} "
        f"{total:6d} "
        f"{status}"
    )